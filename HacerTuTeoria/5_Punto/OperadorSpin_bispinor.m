function Sigma_4 = OperadorSpin_bispinor(index, negativo)
% No sé en que caso puede ser útil este operador...

if ~exist('negativo','var')
      negativo = 1;
end

Sigma = MatrizPauli(index);

ceros = zeros(2);
Sigma_4 = [Sigma, ceros; ceros, negativo * Sigma];