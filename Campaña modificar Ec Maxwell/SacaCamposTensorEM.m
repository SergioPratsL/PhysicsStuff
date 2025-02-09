function [E, B] = SacaCamposTensorEM(T)
% Devuelve E y B a partir del tensor electromagnético
% También vale para el momento dipolar eléctrico y magnético (eso espero)
% Si lo que le pasas es una matriz 4x4 que no es un tensor EM,
% te devolverá una mierda.

E = [T(2,1), T(3,1), T(4,1)];

B = [T(4,3), T(2,4), T(3,2)];