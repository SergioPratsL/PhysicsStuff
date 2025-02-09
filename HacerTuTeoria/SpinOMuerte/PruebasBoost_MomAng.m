
%v = [0.5, 0, 0];  % Fue bien
v = [0.4, -0.3, 0.42];
Lambda = Tensor_boosts(v);

N = [0,0,0];
J = [0,1,0];

N_uv = TensorMomAng(N, J)

%N_uv_boost = Lambda^(-1) * N_uv * Lambda
%N_uv_boost = Lambda * N_uv * Lambda^(-1)
N_uv_boost = Lambda * N_uv * Lambda     % Buena!

[Nf, Jf] = Boost_EM(N, J, v);
N_uv_boost_altern = TensorMomAng(Nf, Jf);

%Ahora cogere vectores que generarían N y J y los transformaré.
X = [0,0,0,1];  % Y
P=[0,-1,0,0]; %Z
N_uv_altern = (X')*P - ((X')*P)';

X_boost = Boost(X, v);
P_boost = Boost(P, v);

% Coincide!
N_uv_boost_altern2 = (X_boost')*P_boost - ((X_boost')*P_boost)'
