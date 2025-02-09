function v_norm = NormalizaCustom(vector, factorNorm)
% Devuelve un vector normalizado pero no a 1 sino lo que tú le digas.

v_norm = vector / norm(vector) * factorNorm;


