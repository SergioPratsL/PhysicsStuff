function [Eo, Bo] = CampoInducido(R, v, q)
% Campo inducido, basado el el potencial de Lienard Wiechert y con unidades
% q es la carga que hace el campoo expresada en Coulombs

% Basado en la expresion de la Wikiedia en ingles cuya formulacion me parece mas elegante:
% http://en.wikipedia.org/wiki/Li%C3%A9nard%E2%80%93Wiechert_potential

% Velocidad de la luz m/s
c = 2.9979 * 10^8;

% Constante de permitividad (Farads por metro)
E_0 = 8.8542 * 10^-12;

R_norm = R / norm(R);

v_norm = v / c;

Sigma = fGamma(v_norm);

Coef_Eo = 1/(4*pi*E_0) * q * (1/Sigma^2) * 1/( 1 - dot(R_norm, v_norm) )^3 * 1 / norm(R)^2;

Eo = Coef_Eo * (R_norm - v_norm);

% Aqui uso v en vez de v_norm para ahorrarme una c que iba fuera
% Bo = Coef_Eo * cross( v, R_norm);  % Comentado 09.10.2015

