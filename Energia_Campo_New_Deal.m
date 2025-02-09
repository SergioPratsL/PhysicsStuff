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
num_iter_V = 30000;     % Toma!

r_unit_V = a0 / 200;    % Precision extrema


% Empezamos por el segundo segmento porque si no se desvirtua todo.
E_pot_acumulada = 0;
E_campo_acumulada = 0;
dens_prob_acum = 0;
iter_actual = 1;
r_ini = 0;
E_prot = 0;
Vol_acum = 0;

while iter_actual <= num_iter_V
    
% Radios centrados en el valor que estamos evaluando
    r_fin = (iter_actual + 0.5) * r_unit_V;
    
    r_centro = iter_actual * r_unit_V;       % Para evaluar

% Volumen del anillo.    
    Vol_anillo = 4 * pi / 3 * ( r_fin^3 - r_ini^3 );
    
    Vol_acum = Vol_acum + Vol_anillo;
    
    expon = exp(-2 * r_centro / a0);
    
% En verdad esta densidad de probabilidad es para el nivel 1_0_0    
    dens_prob_acum_r = (1 - 2 * (r_centro / a0)^2 * expon - 2 * (r_centro / a0) * expon - expon );
    
    E_electron = - cte_electr * dens_prob_acum_r / r_centro^2;
    
    E_proton = cte_electr / r_centro^2;
    
%    E_dif = (Permitividad_vacio / 2) * Vol_anillo * ( (E_proton + E_electron)^2 - E_proton^2 );
%    E_dif =  (Permitividad_vacio / 2) * Vol_anillo * (E_electron)^2; 
    E_dif = (Permitividad_vacio / 2) * Vol_anillo * ( (E_proton + E_electron)^2 - E_proton^2 - E_electron^2);
    
    if iter_actual == 100 || iter_actual == 150 || iter_actual == 2500  % || iter_actual == 1500
        a = 1;
    end
    
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

%dens_prob_total = dens_prob_acum

    
% Control de calidad    
%dens_prob_acum_fin = dens_prob_acum_r 

%Vol_acum_final = Vol_acum

%Vol_acum_esperado = 4 * pi / 3 * ((num_iter_V + 0.5) * r_unit_V)^3