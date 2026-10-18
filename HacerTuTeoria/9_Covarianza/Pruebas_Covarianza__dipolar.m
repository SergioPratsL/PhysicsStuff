% addpath('C:\Fisica\Scripts matlab\HacerTuTeoria\9_Covarianza')

% La densidad de momento magnético y eléctrico se puede obtener mediante
% sencillas matrices, vamos a ver si se transforma de forma covariante.

[sx, sy, sz] = MatricesPauli();

[gt, gx, gy, gz] = MatricesGamma();
[at, ax, ay, az] = MatricesAlfa();
Sx = MatrizSpin_4_4(1);
Sy = MatrizSpin_4_4(2);
Sz = MatrizSpin_4_4(3);
[sx, sy, sz] = MatricesPauli();

spinor_base = [1, 0]; 
v = [0,0,0];
%v = [0.4, 0, 0];
%u = [0, 0, 0];
u = [0, 0.35, 0];


p = fGamma(v) * v;
Phi = DiracSpinorPlainWave(p, spinor_base)

MP = TensorDipolarDirac(Phi.');

%[E_A, B_A] = Boost_EM(E_sric, B_sric, v_ca_rot);
MP_u = BoostTensorEE(MP, u)


% Ahora transformar el spinor y calcular la matriz ahí.
Phi_u = (BoostBispinor(Phi, u));

v_u = Vel_Addition_Law(v,u);
p_u = fGamma(v_u) * v_u;
Phi_u_alt = DiracSpinorPlainWave(p_u, spinor_base);

Phi_u_normalized = Phi_u.' / norm(Phi_u);
Phi_u_alt_normalized = Phi_u_alt / norm(Phi_u_alt);

% Estas dos sí que dan igual, pero no igual que MP_u que es lo que queria
% probar... vaya locur.a
MP_u2 = TensorDipolarDirac(Phi_u)
MP_u2_alt = TensorDipolarDirac(Phi_u_alt.');



