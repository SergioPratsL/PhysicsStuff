% En este script quiero usar la idea que N juega el papel de momento
% orbital "externo" al ser aplicado a Px, Py y Pz y al símbolo L-C
% de manera que N es el momento interno de lo que haya en un punto.

% Lo que quiero averigurar es si el spin siempre será el mismo o no.

clear;

p0 = [0, 0, 0];

% Por definición sólo quiero momento angular intrínseco
N = [0,0,0];   

% El spinor en reposo en Z ya está bien.
spinor = [1, 0];

% Velocidad a la que se mueve el objeto para el laboratorio
%v = [0.4, 0, 0];       % OK
%v = [0, 0, 0.4];       % OK
v = [-0.3, 0.5, 0.4];   % OK

biSpinor = DiracSpinorPlainWave(p0, spinor);

dir_origen = SpinorToVector(biSpinor(1:2));

J = 0.5 * dir_origen;

biSpinor_boost = BoostBispinor(biSpinor, -v)'

spinor_boost_1 = biSpinor_boost(1:2);
spinor_boost_2 = biSpinor_boost(3:4);

dir_spinor1 = SpinorToVector(spinor_boost_1);

dir_spinor2 = SpinorToVector(spinor_boost_2);

dens_onda = norm(spinor_boost_1)^2 + norm(spinor_boost_2)^2;

J_boost = 0.5 * (norm(spinor_boost_1)^2 * dir_spinor1 + norm(spinor_boost_2)^2 * dir_spinor2) / dens_onda;

M_boost = TensorMomAng(N, J_boost);

P = Momento_4(v);

W = Pseudovector_P_L(M_boost, P)



W_ori = [0, J]



W_boost_desde_ori = Boost(W_ori, v)

