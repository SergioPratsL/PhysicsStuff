% addpath('C:\Fisica\Scripts matlab\HacerTuTeoria\9_Covarianza')

% La densidad de momento magnético y eléctrico se puede obtener mediante
% sencillas matrices, vamos a ver si se transforma de forma covariante.

clear

[sx, sy, sz] = MatricesPauli();

[gt, gx, gy, gz] = MatricesGamma();
[at, ax, ay, az] = MatricesAlfa();
Sx = MatrizSpin_4_4(1);
Sy = MatrizSpin_4_4(2);
Sz = MatrizSpin_4_4(3);
[sx, sy, sz] = MatricesPauli();

spinor_base = [1, 0]; 
% OK
% v = [0,0,0];
% u = [0.4, 0, 0];

% OK
% v = [0.4, 0, 0];
% u = [0, 0.4, 0];

% OK
% v = [0, 0, 0.5];
% u = [0, 0.4, 0];

% OK.
 % v = [0.4, 0, 0];
 % u = [0, 0, 0.4];

 % OK.
 % v = [0.4, 0, 0.2];
 % u = [0, 0, 0.4];

 % OK
 % v = [0.4, 0, 0.2];
 % u = [-0.2, 0, 0.4];
 
 % OK - Bug de BoostBispinor resuelto
 % v = [0, -0.4, 0.2];
 % u = [0, 0, 0.01]; 
 
 % OK
 %v = [0, 0.4, 0.2];
 %u = [0, 0, 0.4];
 
 % OK - ¡¡¡VICTORIA DECISIVA!!!
 v = [0.32, -0.27, 0.5];
 u = [0.51, 0.4, -0.1];

% todo ok, pero el Pruebas_Covarianza_Spinors ahora está fallando :(.
% 25.10.2026 19:04. Tanto este ccomo el Pruebas_Covarianza_Spinors.m OK!!
 
p = fGamma(v) * v;
Phi = DiracSpinorPlainWave(p, spinor_base)

% %MP = TensorDipolarDirac(Phi.')
[P, M] = DensidadDipolarDirac(Phi.');
MP_i = TensorEM(P, M)

%MP_u = BoostTensorEE(MP, u)
[P_u, M_u] = Boost_Dipolar(P, M, u);
MP_u = TensorEM(P_u, M_u)

% Ahora transformar el spinor y calcular la matriz ahí, se detecto un
% problema que ya ha sido arreglado.
Phi_u = (BoostBispinor(Phi, u));
%J_4 = CorrienteBispinor(Phi_u.');
% Cambia la corriente si conjugo el vector? Sí, la Y
J_4 = CorrienteBispinor(Phi_u');
J = J_4(2:4)/J_4(1);

v_u = Vel_Addition_Law(v,u);
%p_u = fGamma(v_u) * v_u;
%%Phi_u_alt = DiracSpinorPlainWave(p_u, spinor_base);

%Phi_u_normalized = Phi_u.' / norm(Phi_u);
%Phi_u_alt_normalized = Phi_u_alt / norm(Phi_u_alt);

% Estas dos sí que dan igual, pero no igual que MP_u que es lo que queria
% probar... vaya locura.
%MP_u2 = TensorDipolarDirac(Phi_u)
% ¡Fuera la infamia de conjugar porquee sí!
%[P_u2, M_u2] = DensidadDipolarDirac(conj(Phi_u));
[P_u2, M_u2] = DensidadDipolarDirac(Phi_u);
MP_u2 = TensorEM(P_u2, M_u2)



