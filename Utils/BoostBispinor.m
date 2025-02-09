function BispinorRotado = BoostBispinor(bispinor, v)
% Aplica boosts a bispinors de Dirac

v_norm = v / norm(v);
rapidez = atanh(norm(v));
rapidez_vect = rapidez * v_norm;

expM = MatrizBoostSpinDirac(rapidez_vect);

BispinorRotado = expM * bispinor';