% Este será un script pequeño
clear;

N = [0, 0, 0];
J = [0, 0, 1];


% Prueba 1:  v = [0.4, 0, 0]
%N_boost = N_boost_prec = [0  0.4364  0]

%Prueba 2:   v = [0, 0.4, 0];
%N_boost = [-0.4364  0  0]
%N_boost_prec = [-0.4364  0  0.0004]

%Prueba 3:  v = [0, 0, 0.4];
%N_boost = [0,  0,  0]
%N_boost_prec = 1.0e-03 *[0  -0.4364  0]

v = [0.4, 0, 0];
%v = [0, 0.4, 0];
%v = [0, 0, 0.4];


M_ori = TensorMomAng(N, J)
M_boost = Boost_Tensor_Momento_Angular(N, J, v)
% Si el momento de masa se desplaza, será que la masa se ha desplazado por
% efecto del momento angular.
N_boost = M_boost(1, 2:4)

%Ahora una pequeña precesión.
J_prec = J + [0.001, 0, 0];
M_prec = TensorMomAng(N, J_prec)

M_boost_prec = Boost_Tensor_Momento_Angular(N, J_prec, v)

N_boost_prec = M_boost_prec(1, 2:4)
