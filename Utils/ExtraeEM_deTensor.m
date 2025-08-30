function [E, B] = ExtraeEM_deTensor(T_EM)
% Esta función es una utilidad, ¡No comprueba que el tensor EM sea
% antisimétrico!

E = [T_EM(2, 1), T_EM(3, 1), T_EM(4, 1)];

B = [T_EM(4, 3), T_EM(2, 4), T_EM(3, 2)];
