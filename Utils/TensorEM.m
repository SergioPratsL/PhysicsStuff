function T_EM = TensorEM(E, B)
% Funci�n que monta el tensor electromagnetico a partir de los campos E y B
% Se omiten todas las constantes para trabajar de forma simplificada
% Devuelve el tensor en la forma covariante T^uv 
% Valido tambi�n para los momentos dipolares el�ctrico y magn�tico

T_EM = zeros(4);

T_EM(1, 2) = - E(1);
T_EM(1, 3) = - E(2);
T_EM(1, 4) = - E(3);

T_EM(2, 1) = E(1);
T_EM(3, 1) = E(2);
T_EM(4, 1) = E(3);

T_EM(3, 2) = B(3);
T_EM(4, 2) = -B(2);

T_EM(2, 3) = -B(3);
T_EM(4, 3) = B(1);      % ojo!

T_EM(2, 4) = B(2);      % ojo!
T_EM(3, 4) = -B(1);     % ojo


