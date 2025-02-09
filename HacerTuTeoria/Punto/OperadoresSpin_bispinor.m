function [Sigma_4x, Sigma_4y, Sigma_4z]  = OperadoresSpin_bispinor(negativo)

if ~exist('negativo','var')
      negativo = 1;
end

Sigma_4x = OperadorSpin_bispinor(1, negativo);
Sigma_4y = OperadorSpin_bispinor(2, negativo);
Sigma_4z = OperadorSpin_bispinor(3, negativo);
