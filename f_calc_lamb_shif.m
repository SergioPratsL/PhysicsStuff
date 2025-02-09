% Matlab para calcular el desplazamiento lamb con mi formula

clear;      % Borrado de variables que existieran previamente

% DECLARACION DE CONSTANTES en el Sistema de Referencia Internacional:

% Radio de Bohr:   5.2917721067(12)×10?11 m
a0 = 5.2918 * 10^-11;
a0_por_dos = 2 * a0;

% Carga del electron:  1.602 176 565(35)×10?19  Coulombs
e = 1.6022 * 10^-19;

% Constante de Planck:  6.62606957 × 10?34 J * s
h = 6.6261 * 10^-34;

% Masa del electron: 9.10938356 × 10-31
m_e = 9.10938 * 10^-31;

%Velocidad de la luz: 299792458 m / s
c = 2.9979 * 10^8;
c_cuadr = c^2;

% Constante de Planck 6.62607004 × 10-34 m2 kg / s
c_planck = 6.6261 * 10^-34;
c_planck_reduc = c_planck  / (2*pi);

% Constante que normaliza la onda del nivel 2_0_0 elevada al cuadrado,
% sirve para normalizar la carga.
cte = 1/(32*pi*a0^3);
cte_onda = sqrt( cte );     % La constante anterior 'cte' es la constante de la onda al cuadrado

% Energia del nivel 2: 3.4015 eV, pasarla a J
e_nivel_2 = 3.4015 * e;

Permitividad_vacio = 8.85419 *10^-12;       % Vergonzoso, me lo habia dejado

% Termino energia: es el resultado de hacer - 2 mc^2 * 3.4015e, empleado en
% el propagador
term_ener_cuadr = 2 * e_nivel_2 / (c^2 * m_e);    % Las c^2 a efectos de desglosar todos los factores    

sacar_fft = 1;


% La variable nivel controla si evaluamos la autoenergia para el nivel del
% atomo de hidrogeno '2_0_0' o para el '2_1_0'.
nivel = '2_0_0';

% Numero de coordenadas del mapa de la onda que creare
num_iter_cartesiana = 128;
%num_iter_cartesiana = 16;           % Test
num_iter_cartesiana_mitad = num_iter_cartesiana / 2;

% Segmentos radiales
%r_unit = a0 / 10;
r_unit = a0 * 0.2;        % Test, habria que cambiar los ficheros de potenciales!!
%num_iter_V = 500;
num_iter_V = 300;
r_unit_V = a0 / 10;  % Esto viene del fichero del potencial.


% Segmentos de colatitud
num_iter_ang = 90;
ang_unit = pi / num_iter_ang; 

% Otro error garrafal.
%vol_element = r_unit^3 * (1/num_iter_cartesiana^3);
vol_element_bad = (h /(2*pi))^3 * (1/(r_unit/num_iter_cartesiana))^3;
vol_element = r_unit^3;

% Informar en la consola la fecha de inicio del proceso
fecha_hora = clock;

hora = round(10 * fecha_hora(4) / 10);
minuto = round(10 * fecha_hora(5) / 10);
segundo = round(10 *fecha_hora(6) / 10);

fprintf('Inicio: %d:%d:%d \r\n', hora, minuto, segundo);


% Saber los cosenos para elegir el angulo correcto
% Sacar senos y cosenos para no tener que calcularlos tantas veces
coseno = zeros(1, num_iter_ang);
seno = zeros(1, num_iter_ang);
iter_ang_actual = 1;
while iter_ang_actual <= num_iter_ang
    ang_actual = - pi / 2 + iter_ang_actual * ang_unit;
%    ang_actual =  iter_ang_actual * ang_unit;      % Iteraciones de -pi/2
%    a pi/ 2
    coseno(iter_ang_actual) = cos(ang_actual);    
    seno(iter_ang_actual) = sin(ang_actual);    
    iter_ang_actual = iter_ang_actual + 1;
end



if nivel == '2_0_0'
    fileID = fopen('C:\Datos\Pot_Atm_H_2_0_0_indir.txt','r');
elseif nivel == '2_1_0'
    fileID = fopen('C:\Datos\Pot_Atm_H_2_1_0_indir.txt','r');
else
    fprintf('Nivel de energia erroneo');
    return
end


% FSCANF exige que me cargue las lineas de comentarios que suelo poner
formato = '%f ;'
V = zeros(num_iter_ang, num_iter_V);

cont = 1;
while cont <= num_iter_V
    V_esf = fscanf( fileID, formato);
    V(1:90, cont) = V_esf;
    cont = cont + 1;
end

% Recordar que estos potenciales estan multiplicados por 10^10!!
% Vergonzosamente me habia dejado la permitividad:
% V = V * 10^-10 / (4*pi*Permitividad_vacio);       % Ya no hace falta :).



f_onda = zeros(num_iter_cartesiana, num_iter_cartesiana, num_iter_cartesiana);
f_onda_pot = zeros(num_iter_cartesiana, num_iter_cartesiana, num_iter_cartesiana);

% Variables de chequeo
Onda_cuadr_acum = 0;
Onda_V_cuadr_acum = 0;


% Obtener la funcion de onda y el producto de la funcion de onda por el
% potencial
cont_x = 1;
iter_ang_anterior = 1;      % Var auxiliar.

cont_aux = 1;

l_offset = - (num_iter_cartesiana / 2 + 0.5);

while cont_x <= num_iter_cartesiana
    cont_y = 1;
    x = r_unit * ( cont_x + l_offset );

    while cont_y <= num_iter_cartesiana
        cont_z = 1;
        y = r_unit * ( cont_y + l_offset );
        
        while cont_z <= num_iter_cartesiana;
            
% debug
            if cont_y == 84 && cont_z == 66
                ada = 120;
            end
            
            z = r_unit * ( cont_z + l_offset );       

% Hay magnitudes que son comunes para las dos funciones de onda            
            r = sqrt( x^2 + y^2 + z^2 );
            exp_onda = exp( - r / a0_por_dos );           
            
% Hay que interporal el potencial en funcion de la posicion
            colatitud = z / r;

            encontrado = 0;            
            if colatitud <= seno(1)
                encontrado = 1;
                num_sen_ant = 1; 
                num_sen_post = 1;                      
                iter_ang = 1;
            elseif iter_ang_anterior > 4
                iter_ang = iter_ang_anterior - 3;
            else
                iter_ang = 1;
            end
            
            sen_ant = seno(iter_ang);
            iter_ang = iter_ang + 1;
% Hacemos una especie de Bolzano... bastante cutre pero ahorrara tiempo:
            while encontrado == 0
                sen_act = seno(iter_ang); 
                
                if sen_act >= colatitud && sen_ant <= colatitud
                    encontrado = 1;
                    num_sen_ant = iter_ang - 1; 
                    num_sen_post = iter_ang;                    
                else
                    sen_ant = sen_act;
                    iter_ang = iter_ang + 1;
                end
                
% No hemos conseguido encontrar el angulo, leer desde el principio (deberia pasar pocas veces)                
                if iter_ang > num_iter_ang
                    encontrado = 2;
                    iter_ang = 1;
                end
            end
            
            if encontrado == 2
                
                if colatitud <= seno(1)
                    num_sen_ant = 1;
                    num_sen_post = 1;
                    encontrado = 1;
                elseif colatitud >= seno(num_iter_ang)
                    num_sen_ant = num_iter_ang;
                    num_sen_post = num_iter_ang;
                    encontrado = 1;               
                end
                    
                sen_ant = seno(iter_ang);
                iter_ang = iter_ang + 1;
                while encontrado ~= 1
                    
                    sen_act = seno(iter_ang); 
                
                    if sen_act >= colatitud && sen_ant <= colatitud
                        encontrado = 1;
                        num_sen_ant = iter_ang - 1;
                        num_sen_post = iter_ang;
                    else
                        sen_ant = sen_act;
                        iter_ang = iter_ang + 1;
                    end
                
% Anomalia, forzar pete
                    if iter_ang > num_iter_ang
                        pete = 1 / 0
                    end
                end
            end                
                    
% Determinar el radio, puede hacerse de forma mas directa.
            valor = r / r_unit_V;           % Aplicar los valores del potencial!!
            num_r_ant = floor( valor );
            num_r_post = ceil( valor );
            
            if num_r_post > num_iter_V
                num_r_post = num_iter_V;
            end

            if num_r_ant < 1
                num_r_ant = 1;
            end
            
            
% Guarradilla: para sacar la interpolacion de los angulos me basare en los cosenos y no en los angulos.
            dist_ang = seno(num_sen_post) - seno(num_sen_ant);
            
            if dist_ang == 0    % Excepcion
                V1_r = V(num_sen_post, num_r_ant);
                V2_r = V(num_sen_post, num_r_post);
                                
            else                % Caso general
                V1_r = (seno(num_sen_post) - colatitud) * V(num_sen_post, num_r_ant) + (colatitud - seno(num_sen_ant)) * V(num_sen_ant, num_r_ant);
                V1_r = V1_r / dist_ang;
                
                V2_r = (seno(num_sen_post) - colatitud) * V(num_sen_post, num_r_post) + (colatitud - coseno(num_sen_ant)) * V(num_sen_ant, num_r_post) ;
                V2_r = V2_r / dist_ang;
                
            end
            
% Por fin obtenemos el potencial interpolado!
            r2 = r_unit * num_r_post;
            r1 = r_unit * num_r_ant;

            if r == r2 or r2 == r1
                V_obtenido = V2_r;
            elseif r == r1
                V_obtenido = V1_r;
            else
                V_obtenido = (r2 - r) * V1_r + (r - r1) *  V2_r; 
                V_obtenido = V_obtenido / r_unit;           % Me lo deje, fallo garrafal!
            end
                        
                
            if nivel == '2_0_0'    
            
               f_onda(cont_x, cont_y, cont_z) = cte_onda * (2 - r / a0) * exp_onda;
               
            elseif nivel == '2_1_0'

%               cos_ang = sqrt( (x^2 + y^2) / r^2 );
                cos_ang = colatitud;            % Por usar el rango de -pi/2 a pi/2 (mas facil que cambiar todo lo demas) :).
                f_onda( cont_x, cont_y, cont_z) = cte_onda * (r / a0) * exp_onda * cos_ang;                               
                
            end
        
            f_onda_pot(cont_x, cont_y, cont_z) = f_onda( cont_x, cont_y, cont_z) * V_obtenido;            
            
            Onda_cuadr_acum = Onda_cuadr_acum + f_onda( cont_x, cont_y, cont_z)^2;
            Onda_V_cuadr_acum = Onda_V_cuadr_acum + f_onda_pot( cont_x, cont_y, cont_z)^2;
            
            cont_z = cont_z + 1;
            cont_aux = cont_aux + 1;           
        end
        
        cont_y = cont_y + 1;
    end
    
    cont_x = cont_x + 1;
end


% Matrices f_onda y f_onda_pot terminadas
fecha_hora = clock;

hora = round(10 * fecha_hora(4) / 10);
minuto = round(10 * fecha_hora(5) / 10);
segundo = round(10 *fecha_hora(6) / 10);

fprintf('Matrices f_onda y f_onda_pot terminadas: %d:%d:%d \r\n', hora, minuto, segundo);

return



% Hacemos las FFTs...
f_onda_p = fft( f_onda );

%clear f_onda;

% Conjugar la funcion de onda!!
f_onda_p = conj(f_onda_p);

f_onda_pot_p = fft( f_onda_pot ); 

%clear f_onda_pot;

% Matrices f_onda y f_onda_pot terminadas
fecha_hora = clock;

hora = round(10 * fecha_hora(4) / 10);
minuto = round(10 * fecha_hora(5) / 10);
segundo = round(10 *fecha_hora(6) / 10);

fprintf('FFTs terminadas: %d:%d:%d \r\n', hora, minuto, segundo);


f_unit = (1 / c * m_e) * 1 / (r_unit * num_iter_cartesiana);
f_unit_cuadr = f_unit^2;

% f_unit_cuadr_c_cuadr = f_unit_cuadr * c_cuadr;    % Esto es incorrecto. 


% Variables de chequeo
Onda_p_cuadr_acum = 0;
Onda_V_p_cuadr_acum = 0;


Valor_acumulado = 0;
cont_x = 1;
while cont_x <= num_iter_cartesiana
    cont_y = 1;
    
    while cont_y <= num_iter_cartesiana
        cont_z = 1;    
        
        while cont_z <= num_iter_cartesiana;
      
% Ojo! En la DFT las frecuencias por encima de la mitad son en verdad
% frecuencias negativas!
            if cont_x <= num_iter_cartesiana_mitad
                cont_x_aux = cont_x - 1;                
            else
                cont_x_aux = cont_x - num_iter_cartesiana - 1;   
            end

            if cont_y <= num_iter_cartesiana_mitad
                cont_y_aux = cont_y - 1;                
            else
                cont_y_aux = cont_y - num_iter_cartesiana - 1;   
            end            
            
            if cont_z <= num_iter_cartesiana_mitad
                cont_z_aux = cont_z - 1;                
            else
                cont_z_aux = cont_z - num_iter_cartesiana - 1;   
            end            
            
            
            momentum_cuadr = f_unit_cuadr_c_cuadr * ( (cont_x_aux)^2 + (cont_y_aux)^2 + (cont_z_aux)^2);      
            Propagador = 1 / ( term_ener_cuadr + momentum_cuadr);               
    
%            Propagador = 1;     % Test salvaje (y fallido :P)
            Onda_p_cuadr_acum = Onda_p_cuadr_acum + f_onda_p( cont_x, cont_y, cont_z ) * conj(f_onda_p( cont_x, cont_y, cont_z ));
%            Onda_V_p_cuadr_acum = Onda_V_p_cuadr_acum + f_onda_pot_p( cont_x, cont_y, cont_z ) * conj(f_onda_pot_p( cont_x, cont_y, cont_z ));

                
            Valor = f_onda_p( cont_x, cont_y, cont_z ) * f_onda_pot_p( cont_x, cont_y, cont_z ) * Propagador;
            
            Valor_acumulado = Valor_acumulado + Valor;
        
            cont_z = cont_z + 1;
        end
        
        cont_y = cont_y + 1;        
    end
    
    cont_x = cont_x + 1;
end


% Todo el proceso terminado
fecha_hora = clock;

hora = round(10 * fecha_hora(4) / 10);
minuto = round(10 * fecha_hora(5) / 10);
segundo = round(10 * fecha_hora(6) / 10);

fprintf('Fin del proceso: %d:%d:%d \r\n', hora, minuto, segundo);


% Aplicamos el factor de volumen de cada punto sobre el que hemos integrado
Valor_acumulado = Valor_acumulado * vol_element * e;            % No nos olvidemos de la e!!!

Valor_acumulado = Valor_acumulado / num_iter_cartesiana^2;


% Este c_planck_reduc  se debe a que en los calculos se emplea la energia
% que causaria el autopotencial como fuente del efecto shift pero en verdad
% dicha fuente es la variacion en la funcion de onda en origen que causa el
% auto-potencial
%... Si bien esa variacion no es verdadera... es un efecto oscuro...
Lamb_Shift_Pratiano = Valor_acumulado / c_planck_reduc 

Lamb_Shift_Pratiano_frec = Lamb_Shift_Pratiano / h


Onda_p_cuadr_acum
Onda_cuadr_acum

Ratio = Onda_p_cuadr_acum / Onda_cuadr_acum

%Onda_V_p_cuadr_acum
%Onda_V_cuadr_acum

Norm_Check = Onda_p_cuadr_acum * vol_element_bad

Norm_Check2 = Onda_cuadr_acum * vol_element