% Script que obtiene el potencial para el atomo de hidrogeno en el nivel (2,0,0) 
% integrando por distancias (se pone en duda que la integral que obtuve con
% el metodo directo estuviera bien


% DECLARACION DE CONSTANTES en el Sistema de Referencia Internacional:

% Radio de Bohr:   5.2917721067(12)×10?11 m
a0 = 5.2918 * 10^-11;

% Carga del electron:  ?1.602 176 565(35)×10?19  Coulombs
e = 1.6022 * 10^-19;

% Constante que normaliza la onda del nivel 2_0_0 elevada al cuadrado,
% sirve para normalizar la carga.
cte = 1/(32*pi*a0^3);


fecha_hora = clock;

hora = round(10 * fecha_hora(4) / 10);
minuto = round(10 * fecha_hora(5) / 10);
segundo = round(10 *fecha_hora(6) / 10);

fprintf('Inicio: %d:%d:%d \r\n', hora, minuto, segundo);


fileID = fopen('C:\Datos\Pot_Atm_H_2_0_0_semi_dir.txt','w');
formatSpect_lin = '%6.6f \r\n';
formatSpect = '%6.6f %1s';


% Segmentos radiales
r_unit = a0 / 10;
num_iter_r = 100;


iter_r_actual = 1;


fprintf (fileID, '%22s %d \r\n', 'Numero de iteraciones radiales:', num_iter_r);
fprintf (fileID, '%30s \r\n', 'Valores normalizados por 10^10');


while iter_r_actual <= num_iter_r      % Loop 1: radial ...Test!
    r_actual = r_unit * iter_r_actual;   
    r_int_old = 0;
    iter_rad_int = 1;

    pot_acum = 0;
    
    % Loop para evaluar las contribuciones de cada superficie esferica  
    while iter_rad_int <= num_iter_r      % "int" significa interno
            
        r_int = r_unit * iter_rad_int;
            
        dens_carg = e * cte * exp(- r_int / a0);
        
% Factor de la carga de la esfera: (2 - r/ ao)^2
        factor_rad  = (2 - r_int / a0)^2;        
        
        dens_carg = dens_carg * factor_rad;
        
        Vol = 4 * pi / 3 * (r_int^3 - r_int_old^3);
            
        if r_int <= r_actual    
            pot = dens_carg * Vol / r_actual;
        else
            pot = dens_carg * Vol / r_int;    
        end
        
        pot_acum = pot_acum + pot;
        
        iter_rad_int = iter_rad_int + 1;    
        r_int_old = r_int;
        
    end
    
% Normalizamos el potencial
    pot_norm = pot_acum * 10^10;    
    
    fprintf (fileID, formatSpect_lin, pot_norm);
    
    iter_r_actual = iter_r_actual + 1;
    
end


fclose(fileID);


fecha_hora = clock;

hora = round(10 * fecha_hora(4) / 10);
minuto = round(10 * fecha_hora(5) / 10);
segundo = round(10 *fecha_hora(6) / 10);

fprintf('Fin: %d:%d:%d', hora, minuto, segundo);

            
            