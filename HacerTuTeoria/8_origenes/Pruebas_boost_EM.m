% Pruebas de boosts EM para ver que la fuerza se transforma de forma
% coherente, estas pruebas no incluyen las rotaciones necesarias para que
% las velocidades vistas por todos los SRI sean coherentes!

clear;

v = [0.4, 0, 0];
gamma = fGamma(v);

% En MCU el SRI propio ve una fuerza gamma veces mayor...
%E = [0, 1, 0];
%B = [0, 0, 0];

E = [0.5, 1, -0.5];
B = [2, 1, 0];

% Coinciden
%E = [-0.7, 1, 1.2];
%B = [0.5, 0.6, -0.4];

F = E + cross(v, B)

% La f indica magnitudes en el sistema de la carga, al hacer boost a v
[Ef, Bf] = Boost_EM(E, B, v);

% Porque la carga en su SRI tiene v=0. En realidad 
Ff = Ef

dp = [dot(v,F), F];

dp_f = Boost(dp, v);

% Si dp se obtiene en el lab durante "dt", en el SRI propio es un
% intervalo más corto, luego la fuerza tiene que ser mayor!
Ff_from_dp = dp_f(2:4) * gamma;






