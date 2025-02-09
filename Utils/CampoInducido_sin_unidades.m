function [Eo, Bo] = CampoInducido_sin_unidades(R, v)
% Campo inducido, basado el el potencial de Lienard Wiechert y con unidades
% q es la carga que hace el campoo expresada en Coulombs

% q es 1 en este script simplificado por lo que se quita de los argumentos de la funcion!!

% R es la distancia retardada
% v es la velocidad de la carga

% Basado en la expresion de la Wikiedia en ingles cuya formulacion me parece mas elegante:
% http://en.wikipedia.org/wiki/Li%C3%A9nard%E2%80%93Wiechert_potential

% Velocidad de la luz m/s
% c = 2.9979 * 10^8;
c = 1;  % :P

% Constante de permitividad (Farads por metro)
% E_0 = 8.8542 * 10^-12;
E_0 = 1; % :P

R_norm = R / norm(R);

v_norm = v / c;

Sigma = fGamma(v_norm);

% Coef_Eo = 1/(4*pi*E_0) * q * (1/Sigma^2) * 1/( 1 - dot(R_norm, v_norm) )^3 * 1 / norm(R)^2;
% Fuera constantes!!
Coef_Eo = (1/Sigma^2) * 1/( 1 - dot(R_norm, v_norm) )^3 * 1 / norm(R)^2;

Eo = Coef_Eo * (R_norm - v_norm);

% Aqui uso v en vez de v_norm para ahorrarme una c que iba fuera
% Bo = Coef_Eo * cross( v, R_norm);  % Comentado 09.10.2015

% 09.10.2015. La formula anterior puede que estuviera bien pero prefiero ver claramente la formula de los libros.
Bo = cross( R_norm,  Eo);

