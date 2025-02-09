function [var_Boost, var_Boost_EM] = Boost_vs_BoostEM(Ei, Bi, v_boost)

% Ei y Hi es el campo EM, por tanto dos vectores
% v_boost es el boost que harás para evaluar el campo en otro SRI. Debe
% venir normalizada (c=1)

% Se asume para simplificar que la partícula está quita en SRI origen

% var_Boost es el boost a al 4 vector potencia-fuerza en SRI origen [0, Ei]

% var_Boost_EM es el resultado de evaluar el 4 vector en el SRI destino.

Pot_F_ini = [0, Ei];

var_Boost = Boost( Pot_F_ini, v_boost);

veloc = norm(v_boost);
Sigma = 1 / sqrt(1-veloc^2);

v = -v_boost;

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost);

F_fin = FLorentz(v, Ef, Bf);

var_Boost_EM = [dot(v, F_fin), F_fin] * Sigma; 



