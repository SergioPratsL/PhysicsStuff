function bispinor = DiracSpinorPlainWaveNormalized(p, fi, negativo)
% Obtiene un bispinor de onda plana cuya densidad será gama(v)
% en vez de ser 1 como sugiere la página de la Wikipedia:
% https://en.wikipedia.org/wiki/Dirac_spinor  puesto que una 
% transformación de una onda plana lleva a este por la 
% contracción de longitudes.

% PD: el nombre de la función es un poco confuso :P

% PD2: puede que esto se tenga que ir a la mierda.

if ~exist('negativo','var')
      negativo = 1;
end

bispinor_not_norm = DiracSpinorPlainWave(p, fi, negativo);

gamma = sqrt(1+norm(p)^2);

% El factor de contracción afecta a la densidad de probabilidad, por lo
% tanto a la onda le corresponde su raíz cuadrada.
bispinor = bispinor_not_norm * sqrt(gamma);