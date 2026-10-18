function MatrizBoost = MatrizBoostSpinDirac(rapidez)

factor = rapidez / 2;

[gt, gx, gy, gz] = MatricesGamma();

M = gt * (gx*factor(1) + gy*factor(2) + gz*factor(3));

% 18.10:  MatrizBoost = ExponencialMatriz(M, 0); es la que estaba vigente 
%MatrizBoost = ExponencialMatriz(M, 0);
% ... pero al restaurar la anterior se arregló un fallo de
% C:\Fisica\Scripts matlab\HacerTuTeoria\9_Covarianza\Pruebas_Covarianza_Spinors.m
MatrizBoost = ExponencialMatriz(M, 0).';

end
