clear;

vBA = [0.4, 0, 0];
vBC = [0, 0.4, 0];
%vect_ini_A = [-0.4, 1, 0, 0];
% resultado: vect_ini_x_C = 0    0.9165         0         0

Sigma = fGamma(vBA);
vect_ini_A = [0.4 / Sigma, 0, 1, 0];
% resultado: vect_ini_x_C = 0    0.1600    0.9165         0

vect_B = Boost(vect_ini_A, -vBA);
vect_C = Boost(vect_B, vBC)


vect_ini_C2 = [0.4 / Sigma, 1, 0, 0];
vect_B2 = Boost(vect_ini_C2, -vBC);
vect_A2 = Boost(vect_B2, vBA);


vAC = Vel_Addition_Law(vBC, vBA);
vAC_norm = vAC / norm(vAC);
vCA = Vel_Addition_Law(vBA, vBC);
vCA_norm = vCA / norm(vCA);

angulo = acos(-dot(vAC_norm, vCA_norm)) * 180 / pi

vect_ini_A3 = [norm(vAC), vAC_norm]     % angulo_3 = 0.5041
%vect_ini_A3 = [-norm(vAC), -vAC_norm]  % angulo_3 = 0.5041
% Con vCA quedan residuos que indican que no vamos por buen camino
%vect_ini_A3 = [-norm(vCA), vCA_norm]   
vect_B3 = Boost(vect_ini_A3, -vBA);
vect_C3 = Boost(vect_B3, vBC);

... Son 4 vectores...
angulo_3 = acos(dot(vect_ini_A3(2:4), vect_C3(2:4)) / norm(vect_C3(2:4)) / norm(vect_ini_A3(2:4))) * 180 / pi

vect_ini_A4 = [0, vAC_norm(2), -vAC_norm(1), 0];
vect_B4 = Boost(vect_ini_A4, -vBA);
vect_C4 = Boost(vect_B4, vBC);

angulo_4 = acos(dot(vect_ini_A4(2:4), vect_C4(2:4)) / norm(vect_C4(2:4)) / norm(vect_ini_A4(2:4))) * 180 / pi





