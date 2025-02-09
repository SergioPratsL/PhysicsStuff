% Matlab para calcular el desplazamiento lamb con mi formula

clear;      % Borrado de variables que existieran previamente

% DECLARACION DE CONSTANTES en el Sistema de Referencia Internacional:

% Radio de Bohr:   5.2917721067(12)×10?11 m
a0 = 5.2918 * 10^-11;

% Carga del electron:  1.602 176 565(35)×10?19  Coulombs
e = 1.6022 * 10^-19;
e_cuad = e^2;


Permitividad_vacio = 8.85419 *10^-12; 
cte_electr = e / (4*pi*Permitividad_vacio);

cte_electr_cuad = cte_electr^2;


nivel = '1_0_0'; 

% Constante que normaliza la onda del nivel 2_0_0 elevada al cuadrado,
% sirve para normalizar la carga.
if nivel == '1_0_0';                    
    cte = 1/(pi*a0^3);
else
    cte = 1/(32*pi*a0^3);;
end

  
% Segmentos radiales
num_iter_V = 300;
num_iter_V_mitad = num_iter_V / 2;

r_unit_V = a0 / 10;  % Esto viene del fichero del potencial.
r_unit = r_unit_V;      

% Segmentos de colatitud
num_iter_ang = 90;
ang_unit = pi / num_iter_ang; 


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
    ang_actual = - pi / 2 + (iter_ang_actual - 0.5) * ang_unit;
%    ang_actual =  iter_ang_actual * ang_unit;      % Iteraciones de -pi/2
%    a pi/ 2
    coseno(iter_ang_actual) = cos(ang_actual);    
    seno(iter_ang_actual) = sin(ang_actual);    
    iter_ang_actual = iter_ang_actual + 1;
end


if nivel == '1_0_0'
    fileID = fopen('C:\Datos\Pot_Atm_H_1_0_0_indir.txt','r');
elseif nivel == '2_0_0'
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


% Empezamos por el segundo segmento porque si no se desvirtua todo.
E_pot_acumulada = 0;
E_campo_acumulada = 0;
dens_prob_acum = 0;
iter_actual = 2;
iter_old = 1;
r_ini = 0;
E_prot = 0;

while iter_actual <= num_iter_V
    
% Radios centrados en el valor que estamos evaluando
    r_fin = (iter_actual + 0.5) * r_unit_V;
    
    r_centro = iter_actual * r_unit_V;       % Para evaluar

% Volumen del anillo.    
    Vol_anillo = 4 * pi / 3 * ( r_fin^3 - r_ini^3 );
                    
% Sacamos el campo electrico como el gradiente del potencial.    
% Vital el signo menos!!   
    if iter_actual == 1
        E_electron = - V(45, iter_actual) / r_unit_V;
    else
        E_electron = - (V(45, iter_old) - V(45, iter_actual)) / r_unit_V;
    end
      
    E_proton = cte_electr / r_centro^2;
      
    if iter_actual == 10 || iter_actual == 70 || iter_actual == 150
       iter_actual = iter_actual;
   end
    
% Tratando solo el campo del proton            
    E_dif = (Permitividad_vacio / 2) * Vol_anillo * ( (E_proton + E_electron)^2 - E_proton^2 );
        
    E_campo_acumulada = E_campo_acumulada + E_dif;
            
% Potencial: tener en cuenta no el infame autopotencial sino el potencial del electron y tener en cuenta la densidad de probabilidad y el volumen
    if nivel == '1_0_0'; 
        dens_prob = cte * exp( - 2 * r_centro / a0);
    elseif nivel == '2_0_0'
        dens_prob = cte * ( 2 - r_centro / a0 )^2 * exp( - r_centro / a0);
    end
    
    E_pot = - cte_electr / r_centro * e * Vol_anillo * dens_prob;
    
    E_pot_acumulada = E_pot_acumulada + E_pot;

    dens_prob_acum = dens_prob_acum + Vol_anillo * dens_prob;
    
    E_prot = E_prot + E_electron * e * r_unit_V;
    
    iter_old = iter_actual;
    r_ini = r_fin;
    iter_actual = iter_actual + 1;
end
            
E_pot_total = E_pot_acumulada
E_pot_total_eV = E_pot_total / e

E_campo_total = E_campo_acumulada
E_campo_total_eV = E_campo_total / e

E_prot_total = E_prot
E_prot_total_eV = E_prot_total / e

dens_prob_total = dens_prob_acum