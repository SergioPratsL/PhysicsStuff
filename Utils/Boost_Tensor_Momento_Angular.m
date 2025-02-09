function M_uv = Boost_Tensor_Momento_Angular(N, J, v)
% Transforma el momento tensor de momento angular a partir de sus
% componentes base, que son N el momento de masa y L el tensor de momento
% angular.

[Nf, Jf] = Boost_EM(N, J, v);

M_uv = TensorMomAng(Nf, Jf);