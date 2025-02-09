function [var_Boost, var_Boost_EM, var_ori] = Boost_vs_BoostEM_FPrats(v_ini, Ei, Bi, v_boost)
% Este matlab solo evalua el trabajo hecho por el nuevo término

% Ei y Hi es el campo EM, por tanto dos vectores
% v_boost es el boost que harás para evaluar el campo en otro SRI. Debe
% venir normalizada (c=1)

% Al contrario que en el anterior, en este caso no tiene que estar la
% partícula quieta en origen, en verdad que esté quita en el origen me
% joderia vivo...

% var_Boost es el boost a al 4 vector potencia-fuerza en SRI origen [0, Ei]

% var_Boost_EM es el resultado de evaluar el 4 vector en el SRI destino.

% Deben incluirse los terminos de Lorentz o no funcionara! cte = 1!!

F_ini = FPrats(v_ini, Ei, Bi);
F_ini_aux = FLorentz(v_ini, Ei, Bi);

F_ini = F_ini + F_ini_aux;

Pot_F_ini = [dot(F_ini, v_ini), F_ini(1), F_ini(2), F_ini(3)];

var_ori = Pot_F_ini;        % Para testing. 

Dilat_ini = 1 / sqrt(1-norm(v_ini)^2);

var_Boost = Boost( Pot_F_ini, v_boost);

v_fin = Vel_Addition_Law(v_ini, v_boost);

veloc = norm(v_fin);
Dilat_fin = 1 / sqrt(1-veloc^2);


[Ef, Bf] = Boost_EM(Ei, Bi, v_boost);

F_fin = FPrats(v_fin, Ef, Bf);
F_fin_aux = FLorentz(v_fin, Ef, Bf);

F_fin = F_fin + F_fin_aux;

var_Boost_EM = [dot(v_fin, F_fin), F_fin] * Dilat_fin / Dilat_ini; 



