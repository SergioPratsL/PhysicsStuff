% Script que obtiene el potencial para el atomo de hidrogeno en el nivel (2,0,0) 
% integrando por distancias y angulos solidos (a la usanza del algoritmo 
% necesario para el nivel (2, 1, 0)

clear

% DECLARACION DE CONSTANTES en el Sistema de Referencia Internacional:

% Radio de Bohr:   5.2917721067(12)×10?11 m
a0 = 5.2918 * 10^-11;

% Carga del electron:  ?1.602 176 565(35)×10?19  Coulombs
e = 1.6022 * 10^-19;

h = 6.6261 * 10^-34;

Permitividad_vacio = 8.85419 *10^-12;       % Vergonzoso, me lo habia dejado


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
    coseno(iter_ang_actual) = cos(ang_actual);
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
V = zeros(num_iter_ang, num_iter_r);

cont = 1;
while cont <= num_iter_r
    V_esf = fscanf( fileID, formato);
    V(1:90, cont) = V_esf;
    cont = cont + 1;
end

% Recordar que estos potenciales estan multiplicados por 10^10!!
% Vergonzosamente me habia dejado la permitividad:
V = V * 10^-10 / (4*pi*Permitividad_vacio);


Carga_tot = 0;
Carga_tot_aux = 0;  % Test
r_old = 0;
r_old_vol = 0;
Energia_total = 0;
iter_r_actual = 1;
Vol_tot_tot = 0;

% Evaluar todos los puntos requiere dos loops y cada punto requiere otros
% dos.

while iter_r_actual <= num_iter_r      % Loop 1: radial 
    iter_ang_actual = 1;
    r_actual = r_unit * iter_r_actual;   
    r_actual_vol = r_actual + r_unit / 2;
    

    Vol_tot = 0;
    Fact_tot = 0;
    if nivel == '2_1_0'
    
        dens_carg_radial = e * cte * (r_actual / a0)^2 * exp(- r_actual / a0);
        
        seno_old = -1;
        while iter_ang_actual <= num_iter_ang   % Loop 2: latitudinal 

            ang_actual = - pi / 2 + iter_ang_actual * ang_unit;
        
% Obtener el volumen del anillo        
            if iter_ang_actual == 1
                dif_alt = seno(1) + 1;
            else
                iter_aux = iter_ang_actual - 1;
                dif_alt = seno(iter_ang_actual) - seno(iter_aux);
            end

            longitud = coseno(iter_ang_actual);         % Hay que sacar el 2pi 
        
% En el limite seria integrar el coseno al cuadrado.  
%           Vol = 4*pi/3 * dif_alt * longitud * (r_actual^3 - r_old_vol^3) / 2;
%            Vol = 4 * pi/3 * coseno(iter_ang_actual)^2 * (r_actual^3 - r_old_vol^3) * 2 / num_iter_ang;

            AngSolid = 2 * pi * (seno(iter_ang_actual) - seno_old);
            Vol = AngSolid * (r_actual_vol^3 - r_old_vol^3) / 3;
            
% Obtener la densidad de carga        
            dens_carg = dens_carg_radial * seno(iter_ang_actual)^2;           
        
            
            Energia_total = Energia_total + Vol * dens_carg * V(iter_ang_actual, iter_r_actual);       
            Carga_tot = Carga_tot + dens_carg * Vol;
        
            seno_old = seno(iter_ang_actual);
            iter_ang_actual = iter_ang_actual + 1;
            
            Vol_tot = Vol_tot + Vol;
            Fact_tot = Fact_tot + Vol * seno(iter_ang_actual-1)^2;
            
        end

% Necesario para arreglar mis limitaciones geometricas...        
        Ratio = Vol_tot / (4 * pi/3 * (r_actual_vol^3 - r_old_vol^3));
        Ratio2 = Fact_tot / Vol_tot;
        Carga_tot_aux = Carga_tot_aux + (4 * pi/9 * (r_actual_vol^3 - r_old_vol^3)) * dens_carg_radial;
%        Vol_tot_tot = Vol_tot_tot + Vol_tot;
         Vol_tot_tot = Vol_tot_tot + (4 * pi/3 * (r_actual_vol^3 - r_old_vol^3));

        if Ratio2 < 0.33
            Ratio2 = Ratio2;
        end

% 2_0_0        
    else
        
        dens_carg = e * cte * exp(- r_actual / a0) * (2 - r_actual / a0)^2;
        Vol = 4*pi/3 * (r_actual_vol^3 - r_old_vol^3);
    
        Carga_tot = Carga_tot + Vol * dens_carg;
        Energia_total = Energia_total + Vol * dens_carg * V(45, iter_r_actual);       
        
    end
    
    r_old = r_actual;
    r_old_vol = r_actual_vol;
    iter_r_actual = iter_r_actual + 1;
end



fclose(fileID);


fprintf ('Terminado');

Energia_total_e = Energia_total / e

Carga_tot_e = Carga_tot / e

Carga_tot_aux_e = Carga_tot_aux / e

Vol_teorico = 4*pi/3 * ( r_unit * num_iter_r )^3;

Ratio_Vol = Vol_tot_tot / Vol_teorico