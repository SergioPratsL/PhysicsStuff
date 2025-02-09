function M_4 = Momento_4(v);
% Da el 4-momento asociado a una partícula que se mueva con velocidad v

Sigma = fGamma(v);

M_4 = Sigma * [1, v];

