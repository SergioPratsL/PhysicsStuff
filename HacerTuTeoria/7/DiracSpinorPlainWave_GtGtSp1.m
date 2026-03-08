function DSp = DiracSpinorPlainWave_GtGtSp1(p, Sp1)
%Como DiracSpinorPlainWave_GtGt pero en vez de passar fi, pasamos el Sp1 y
%con él podemos deducir fi y cascar al Sp2 (E-1)*fi
% El 1 de "E-1" representa la masa.
% el DiracSpinorPlainWave_GtGtSp2 fue una versión anterior, idéntica pero
% el spin fuerte era el 2, y no hace falta cambiar y hacerlo menos
% intuitivo.

% https://en.wikipedia.org/wiki/Dirac_spinor

if (norm(p) == 0)
    DSp = [Sp1, [0,0]];
    return
end

norm_p = norm(p);

P_Dot_PauliVector_norm = (PauliVectorEscalarProd(p / norm_p))';

Sp2 = P_Dot_PauliVector_norm * Sp1.';

E = sqrt(norm_p^2 + 1);

%p_dot_fi = P_Dot_PauliVector * fi.';

factor_Sp2 = (E-1)/ norm_p;

DSp = [Sp1.'; factor_Sp2 * Sp2].';

DSp = DSp / norm(DSp);



