function DSp = DiracSpinorPlainWave_GtGt(p, fi)
% Esto da las ondas planas asumiendo que a la masa también le pertoca un
% término gt de forma que tiene la misma matriz que la energía y el
% potencial EM.

% https://en.wikipedia.org/wiki/Dirac_spinor

% p sería el momento de la onda y phi un two-spinor que apuntará en cierta
% dirección.

if (norm(p) == 0)
    DSp = [[0,0], fi];
    return
end

P_Dot_PauliVector = PauliVectorEscalarProd(p);

E = sqrt(norm(p)^2 + 1);

p_dot_fi = P_Dot_PauliVector * fi.';

DSp = [(E-1)*fi.'; p_dot_fi].';

DSp = DSp / norm(DSp);



