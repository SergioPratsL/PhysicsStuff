clear

% DECLARACION DE CONSTANTES en el Sistema de Referencia Internacional:

% Radio de Bohr:   5.2917721067(12)×10?11 m
a0 = 5.2918 * 10^-11;

% Carga del electron:  ?1.602 176 565(35)×10?19  Coulombs
e = 1.6022 * 10^-19;

h = 6.6261 * 10^-34;

Permitividad_vacio = 8.85419 *10^-12;       % Vergonzoso, me lo habia dejado


cte_V = e / (Permitividad_vacio * 4 * pi );

nivel = '2_1_0';


% Constante que normaliza la onda del nivel 2_0_0 elevada al cuadrado,
% sirve para normalizar la carga.
if nivel == '2_0_0';
    cte = 1/(32*pi*a0^3);
else
    cte = 1/(32*pi*a0^3);
end

% Segmentos radiales
r_unit = a0 / 10;
%num_iter_r = 120;
num_iter_r = 500;



Carga_tot = 0;
r_old = 0;
r_old_vol = 0;
Energia_total = 0;
iter_r_actual = 1;

% Evaluar todos los puntos requiere dos loops y cada punto requiere otros
% dos.

while iter_r_actual <= num_iter_r      % Loop 1: radial 
    iter_ang_actual = 1;
    r_actual = r_unit * iter_r_actual;   
    
    r_actual_vol = r_actual + r_unit/2;
   
    
    if nivel == '2_0_0'
        dens_carg = e * cte * exp(- r_actual / a0) * (2 - r_actual / a0)^2;
    else
% La division por 6 es el resultado de integrar el factor sen^2(Rho).        
        dens_carg = e * cte * exp(- r_actual / a0) * (r_actual / a0)^2 / 6;
    end    
    
    Vol_aux = 4*pi/3 * (r_actual_vol^3 - r_old_vol^3);
    dens_carga_aux = exp(- r_actual / a0) * (2 - r_actual / a0)^2;     % Solo para el 2_0_0
    
    Energia_total = Energia_total + Vol_aux * dens_carg * cte_V / r_actual;
    Carga_tot = Carga_tot + Vol_aux * dens_carg;
    
    r_old = r_actual;
    r_old_vol = r_actual_vol;
    iter_r_actual = iter_r_actual + 1;
end




Energia_total_e = Energia_total / e

Carga_tot_e = Carga_tot / e