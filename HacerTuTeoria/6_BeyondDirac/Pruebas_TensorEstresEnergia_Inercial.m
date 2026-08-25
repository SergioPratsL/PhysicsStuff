% La idea es hacer boosts del tensor de estress energia, y ver cómo quedan

m = 1;      % Siempre

v = [0,0,0];    % Siempre (hasta aquí llegamos.

%u = [0.4, 0, 0];
u = [0.4, -0.5, 0.28];      % Parece que esto va bien (con v = [0,0,0])


EE = TensorEstresEnergia_Inercial(m, v);

EE_boost = BoostTensorEE(EE, u)

% Esto es totalmente fullero, pero necesario, y requiere por tanto saber el
% SRI original, lo bueno es que en el camino directo no lo requiere, es el
% camino que debe prevalecer (es la contracción del eje u):
term_contrac = fGamma(u);
EE_boost_2 = TensorEstresEnergia_Inercial(m*term_contrac, -u);


% Extender esto al tensor E-E EM para ver que el BoostTensorEE(EE, u) va
% bien
E = [11, 13, 17];
B = [3, 5, 11];

F = TensorEM(E, B);
T = Tensor_Estres_Energia_EM(F);

u = [0.4, 0.4, -0.4];
T_boost_directo = BoostTensorEE(T, u)

[E_boost, B_boost] = Boost_EM(E, B,u);

F_boost = TensorEM(E_boost, B_boost);
T_boost_indirecto = Tensor_Estres_Energia_EM(F_boost)
