function DSp = DiracSpinorPlainWave(p, fi, negativo)
% DSp significa Dirac spinor, pero también es un homenage a la legendaria
% aplicación
% Masa obviamente normalizada a 1.

% https://en.wikipedia.org/wiki/Dirac_spinor

% p sería el momento de la onda y phi un two-spinor que apuntará en cierta
% dirección.

if ~exist('negativo','var')
      negativo = 1;
end

 
P_Dot_PauliVector = PauliVectorEscalarProd(p);

E = sqrt(norm(p)^2 + 1);

factor_norm = sqrt((E+1)/2);

SmallBispinor = ( P_Dot_PauliVector * fi.' ) / (E+1);

if negativo ~= -1
    DSp = factor_norm * [fi, SmallBispinor.'];
else
    DSp = factor_norm * [SmallBispinor' ,fi];
end




