% Package: prueba geodesica nivel 1

% El parametro a tratar es la velocidad y la posicion de B
% A esta quito y en el centro

% B esta sobre el plano XY con velocidad en direccion X pero
% por coherencia puedo usar vectores

% R es la distancia de B a A
% v es la velocidad de B


% Prueba 1
R = [0,1,0];
v = [0.4, 0, 0];





[Dif_Pot_A, Dif_Pot_B, Rad] = VarPotGeodesico_lv1(R, v)


Dif_Tot = Dif_Pot_A + Dif_Pot_B + Rad

