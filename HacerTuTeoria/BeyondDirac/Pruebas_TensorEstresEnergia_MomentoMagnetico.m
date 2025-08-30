% La idea es hacer boosts del tensor de estress energia, y ver cómo quedan


%u = [0.4, 0, 0];            % OK
u = [0.4, -0.5, 0.28];     %

E = [11, 13, 17];
B = [2, 3, 5];

F = TensorEM(E,B);

% Densidades de moimento
p = [0, 0.1, 0];
m = [0, 0, 1];

% El tensor de densidades dipolares es como el tensor EM
M = TensorEM(p, m);

EE = Tensor_Estres_Energia_Dipolar(F, M)

EE_boost = BoostTensorEE(EE, u)


[E_boost, B_boost] = Boost_EM(E, B,u);
[p_boost, m_boost] = Boost_EM(p, m,u);

F_boost = TensorEM(E_boost, B_boost);
M_boost = TensorEM(p_boost, m_boost);

EE_boost_indirecto = Tensor_Estres_Energia_Dipolar(F_boost, M_boost)

dif_EE = EE_boost - EE_boost_indirecto

