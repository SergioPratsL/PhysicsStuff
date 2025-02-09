function MatrizBoost = MatrizBoostSpinDirac(rapidez)

factor = rapidez / 2;

[gt, gx, gy, gz] = MatricesGamma();

M = gt * (gx*factor(1) + gy*factor(2) + gz*factor(3));

MatrizBoost = ExponencialMatriz(M, 0);
%MatrizBoost = ExponencialMatriz(M, 0).';

end
