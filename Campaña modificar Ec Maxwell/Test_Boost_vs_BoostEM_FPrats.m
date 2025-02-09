
% Caso 1                       --> OK
%v_ini = [0.4, 0, 0];
%Ei = [10,0,0];
%Hi = [0,0,0];
%v_boost = [0.2, 0, 0];


% Caso 2  (añado Ey y Ez)       --> OK
%v_ini = [0.4, 0, 0];
%Ei = [10,4,5];
%Hi = [0,0,0];
%v_boost = [-0.2, 0, 0];


% Caso 3  (añado B)       --> OK
%v_ini = [0.4, 0, 0];
%Ei = [12,-7,5.5];
%Hi = [-1, 1.4, 2];
%v_boost = [-0.2, 0, 0];

% Caso 3B  (todo en Y en vez de X)     --> OK
%v_ini = [0.0, -0.2, 0];
%Ei = [12,-7,5.5];
%Hi = [-1, 1.4, 2];
%v_boost = [0, 0.3, 0];


% Caso 4  (boost en otra direccion)       --> KO
v_ini = [0.4, 0, 0];
Ei = [12,-7,5.5];
%Hi = [-1, 1.4, 2];
Hi = [0, 0, 0];
%%v_boost = [0, 0.3, 0];
v_boost = [0, 0.01, 0];     % miniboost




[var_Boost, var_Boost_EM, var_ori] = Boost_vs_BoostEM_FPrats(v_ini, Ei, Hi, v_boost);

dif_Boost = var_Boost - var_ori

dif_Boost_EM = var_Boost_EM - var_ori




