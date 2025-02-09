% Script que obtiene el potencial para el atomo de hidrogeno en el nivel (2,0,0) 
% integrando por distancias y angulos solidos (a la usanza del algoritmo 
% necesario para el nivel (2, 1, 0)


% DECLARACION DE CONSTANTES en el Sistema de Referencia Internacional:

% Radio de Bohr:   5.2917721067(12)×10?11 m
a0 = 5.2918 * 10^-11;

% Carga del electron:  ?1.602 176 565(35)×10?19  Coulombs
e = 1.6022 * 10^-19;

% Constante que normaliza la onda del nivel 1_0_0 elevada al cuadrado,
% sirve para normalizar la carga.
cte = 1/(pi * a0^3);

Permitividad_vacio = 8.85419 *10^-12;       % Vergonzoso, me lo habia dejado


fileID = fopen('C:\Datos\Pot_Atm_H_1_0_0_indir.txt','w');
formatSpect_lin = '%6.6f \r\n';
formatSpect = '%6.6f %1s';


% Segmentos radiales
r_unit = a0 / 10;
%num_iter_r = 120;
%num_iter_r = 500;
num_iter_r = 300;

% Segmentos de colatitud
num_iter_ang = 90;
ang_unit = pi / num_iter_ang; 


% Sacar senos y cosenos para no tener que calcularlos tantas veces
seno = zeros(1, num_iter_ang);
coseno_cuad = zeros(1, num_iter_ang);

iter_ang_actual = 1;
while iter_ang_actual <= num_iter_ang

    ang_actual = - pi / 2 + (iter_ang_actual - 0.5) * ang_unit;
    seno(iter_ang_actual) = sin(ang_actual);
    cos_cuad(iter_ang_actual) = cos(ang_actual)^2;
    
    iter_ang_actual = iter_ang_actual + 1;
end


fecha_hora = clock;

hora = round(10 * fecha_hora(4) / 10);
minuto = round(10 * fecha_hora(5) / 10);
segundo = round(10 *fecha_hora(6) / 10);

fprintf('Inicio: %d:%d:%d \r\n', hora, minuto, segundo);

fprintf (fileID, '%22s %d \r\n', 'Numero de iteraciones radiales:', num_iter_r);
fprintf (fileID, '%28s %d \r\n', 'Numero de iteraciones por colatitud', num_iter_ang);
fprintf (fileID, '%30s \r\n', 'Valores normalizados por 10^10');

iter_r_actual = 1;

% Evaluar todos los puntos requiere dos loops y cada punto requiere otros
% dos.

while iter_r_actual <= num_iter_r      % Loop 1: radial 
    iter_ang_actual = 1;
    r_actual = r_unit * iter_r_actual;   
    
    r_actual_cuad = r_actual^2;
    
    while iter_ang_actual <= num_iter_ang   % Loop 2: latitudinal 
        iter_rad_int = 1;
        r_int_old = 0;      % Variable auxiliar para el tercer loop
        ang_actual = - pi / 2 + iter_ang_actual * ang_unit;
        
        pot_integrado_tot = 0;
        r_int_vol_old = 0;
        
        
% Loop para evaluar las contribuciones de cada superficie esferica  
        while iter_rad_int <= num_iter_r      % Loop 3: contribuciones de cada radio, "int" significa interno
            iter_ang_int = 1;

% Variable para calcular el potencial total            
            pot_integrado_radio = 0;
            
% El valor de la carga es a nivel de superficie esferica
            r_int = r_unit * iter_rad_int;
            r_int_vol = r_unit * (iter_rad_int + 0.5);
            
            dens_carg = e * cte * exp(- 2* r_int / a0);
% Volumen que corresponde a esta "superficie esferica"
            Vol_rad = (r_int_vol^3 - r_int_vol_old^3) / 3;
% NOTA: se usa 2pi/3 en vez de 4pi/3 porque luego en el angulo solido
% aparece un factor 2.    

            r_int_cuad = r_int^2;       % Para mejorar la eficacia

% Loop para evaluar las contribuciones de cada colatitud, pero colatitud rotada para que el punto que estemos tratando sea el "eje z"   
            while iter_ang_int <= num_iter_ang  % Loop 4: contribuciones de cada colatitud "inclinada" sobre el eje "r"
                ang_int = - pi / 2 + ang_unit * iter_ang_int;
                
% Distancia de cada anillo, elemento critico
                dist = sqrt( r_int_cuad + r_actual_cuad - 2 * r_int * r_actual * seno(iter_ang_int));
                
% Cuando evaluas un punto sobre si mismo da infinito, pero sabemos que esa contribucion es despreciable por lo que podemos descartarla                
                if dist == 0
                    dist = r_int;       % Tosca manera de descartar este valor
                end
                
% Superficie que corresponde a cada 'aro'    
                if iter_ang_int == 1
                    AngSolid = 2 * pi * (seno(1) + 1);
                else
                    AngSolid = 2 * pi * (seno(iter_ang_int) - seno(iter_ang_int - 1));
                end

% Da 2 al final de cada loop :)                
%                superf_acum = superf_acum + superf;            
                
% Carga de la linea: al ser el nivel de energia 2_0_0 no es necesario obtenerla!!!                
                
                pot_integrado_radio = pot_integrado_radio + AngSolid / dist;
                
                iter_ang_int = iter_ang_int + 1;    
            end
    
% Aplicamos aqui al potencial los factores de volumen y densidad de carga "radial", tras el loop para hacerlo solo una vez :)
            pot_integrado_radio = pot_integrado_radio * Vol_rad * dens_carg;
            
            pot_integrado_tot = pot_integrado_tot + pot_integrado_radio;
            
% Aplicamos el factor de normalizacion!
            
            iter_rad_int = iter_rad_int + 1;    
            r_int_old = r_int;
            r_int_vol_old = r_int_vol;
        end
        
% Normalizamos el potencial
        pot_integrado_tot  = pot_integrado_tot / (4*pi*Permitividad_vacio);
        
        if iter_ang_actual == num_iter_ang
            fprintf (fileID, formatSpect_lin, pot_integrado_tot);
        else
            fprintf (fileID, formatSpect, pot_integrado_tot, ';');
        end
        
        iter_ang_actual = iter_ang_actual + 1;
    end
    
    iter_r_actual = iter_r_actual + 1;
end



fclose(fileID);


fecha_hora = clock;

hora = round(10 * fecha_hora(4) / 10);
minuto = round(10 * fecha_hora(5) / 10);
segundo = round(10 * fecha_hora(6) / 10);

fprintf('Fin: %d:%d:%d', hora, minuto, segundo);

fprintf ('Terminado');