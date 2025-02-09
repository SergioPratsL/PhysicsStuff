function spin_boosted = BoostSpin( spin, v);

v_norm = v / norm(v);
xi = atanh(norm(v));
xi_vect = xi * v_norm;

% Se usa el -1i para aprovechar la función MatrizRotacionSpin.
% expM = MatrizRotacionSpin(xi_vect, -1i);

expM = MatrizBoostSpinDirac(xi_vect, -1i);

spin_boosted = expM * spin';