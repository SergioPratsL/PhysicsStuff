% Script que obtiene el potencial para el atomo de hidrogeno en el nivel (2,0,0) 
% integrando por distancias

% DECLARACION DE CONSTANTES en el Sistema de Referencia Internacional:

% Radio de Bohr:   5.2917721067(12)×10?11 m
a0 = 5.2918 * 10^-11;

% Carga del electron:  ?1.602 176 565(35)×10?19  Coulombs
e = 1.6022 * 10^-19;

% Constante que normaliza la onda del nivel 2_0_0 elevada al cuadrado,
% sirve para normalizar la carga.
cte = 1/(32*pi*a0^3);



fileID = fopen('C:\Datos\Pot_Atm_H_2_0_0_dir.txt','w');
formatSpect = '%6.6f \r\n';

% Segmentos
r_unit = a0 / 10;
num_iter = 100;

iter_actual = 1;

fprintf (fileID, '%22s %d \r\n', 'Numero de iteraciones:', num_iter);
fprintf (fileID, '%30s \r\n', 'Valores normalizados por 10^10');


while iter_actual <= num_iter
    
    r_actual = r_unit * iter_actual;
    
% Carga interna: 
    var1_int = 8 * a0^3;            % Integral evaluada en el 0
    var2_int = exp(- r_actual / a0) * (8 * a0^3 + 8 * a0^2 * r_actual + 4 * a0 * r_actual^2 + r_actual^4 / a0);
    
    Q_INT = e * cte * (var1_int - var2_int);
    
    V_INT = Q_INT / r_actual;
    
    V_EXT = e * cte * exp(- r_actual / a0) * (2 * a0^2 + 2 * a0 * r_actual - r_actual^2 + r_actual^3 / a0);
    
    V = V_INT + V_EXT;
    
    V_norm = V * 10^10;
    
    
    fprintf (fileID, formatSpect, V_norm);
    iter_actual = iter_actual + 1;
    
%     if iter_actual > 100
%         fprintf ('%15s \n', 'Alterta!!');
%     end
    
end

fclose(fileID);