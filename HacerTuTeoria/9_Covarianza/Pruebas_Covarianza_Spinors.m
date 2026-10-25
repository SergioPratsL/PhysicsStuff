% addpath('C:\Fisica\Scripts matlab\HacerTuTeoria\9_Covarianza')

% La idea es simplemente transformar un spinor de A --> B --> C
clear;

[sx, sy, sz] = MatricesPauli();

[gt, gx, gy, gz] = MatricesGamma();
[at, ax, ay, az] = MatricesAlfa();
Sx = MatrizSpin_4_4(1);
Sy = MatrizSpin_4_4(2);
Sz = MatrizSpin_4_4(3);
[sx, sy, sz] = MatricesPauli();

spinor_base = [1, 0]; 
%OK
% Velocidad de la onda vista desde A
% v_i = [0,0,0];
% v_ba = [0.4, 0, 0];
% v_bc = [0, -0.4, 0];

% OK
% v_i = [0,0,0];
% v_ba = [0.4, 0.1, 0];
% v_bc = [0.2, -0.36, 0];

% OK
% v_i = [0.25,0,0];
% v_ba = [-0.4, 0.1, 0];
% v_bc = [0.2, -0.46, 0];

% OK
% v_i = [0,0,0];
% v_ba = [0.4, 0, 0];
% v_bc = [0, 0, 0.4];

% OK
% v_i = [0,0,0];
% v_ba = [0, 0, 0.5];
% v_bc = [0.4, 0, 0];

% FALLA ¡¿?! --> OK...
% v_i = [0,0,0];
% v_ba = [0, 0, 0.4];
% v_bc = [0, 0.4, 0];

% KO --> Fue bien!!
% v_i = [0,0,0];
% v_ba = [0.4, 0.1, 0.5];
% v_bc = [0.1, -0.29, 0];

% OK --> Victoria!!!
v_i = [0.24, -0.32, 0.2];
v_ba = [0.4, 0.1, 0.5];
v_bc = [0.1, -0.29, 0];

v_ca = Vel_Addition_Law(v_ba, v_bc);
v_ac = Vel_Addition_Law(v_bc, v_ba);

% Para controlar que sé lo que hago.
%... Si me costó con esto, lógico que me cueste más con spinors
v_test = [1, 1.2, 0.9, -0.6];
%v_test = [1, 0,0,0]
%v_test = [0, 1.2, 0.9, -0.6]

v_test_b = Boost(v_test, -v_ba);
v_test_c = Boost(v_test_b, v_bc);
v_test_rot = Boost(v_test_c, v_ca);

RW = RotacionWigner(v_ba, v_bc);
[dir_rot, angulo] = ExtraeAnguloDeMatriz(RW);

cos_v = dot(v_ca, -v_ac) / norm(v_ac) / norm(v_ca);
angulo_2 = acos(cos_v);
dif_ang = angulo - angulo_2;

v_test_again = (RW * v_test_rot.').';

dif_v_test = v_test - v_test_again

% Pruebas del diario 5 que no se llegaron a hacer:
p0 = [0, 0, 0];
 
% % % Parece que la rotación en la Y no es diferente que en la X
% % spin_base = [1, 0];   
% % bispinor = DiracSpinorPlainWave(p0, spin_base)
% % eje = [1, 1, 0] / sqrt(2);
% % radianes = pi/2;
% % bispinor_rot = RotacionBispinor(bispinor, eje, radianes)
% % sp_rot = bispinor_rot(1:2);
% % dir_rotada = SpinorToVector(sp_rot')
% 
% spin_base = [1, 1] / sqrt(2);   
% dir_ini =  SpinorToVector(spin_base)
% bispinor = DiracSpinorPlainWave(p0, spin_base);
% eje = [0, 1, 1] / sqrt(2);
% radianes = pi/2;
% bispinor_rot = RotacionBispinor(bispinor, eje, radianes);
% sp_rot = bispinor_rot(1:2);
% dir_rotada = SpinorToVector(sp_rot')
% 
% % ~Pruebas del diario 5 que no se llegaron a hacer:
% % ... la verdad es que mucho más debería hacer, estas pruebas
% % solo sirven para bichus en reposo...

p_i = fGamma(v_i) * v_i;
Phi_a = DiracSpinorPlainWave(p_i, spinor_base)

% Ojo, si no conjugas eel Phi_b, o sea, si haces Phi_b.', da mal!
Phi_b = BoostBispinor(Phi_a, -v_ba);
Phi_c = BoostBispinor(Phi_b.', v_bc);
Phi_a_again_rot = (BoostBispinor(Phi_c.', v_ca)).';

Phi_a_again = (RotacionBispinor(Phi_a_again_rot, dir_rot, -angulo)).'

n1 = norm(Phi_a);
n2 = norm(Phi_a_again);



