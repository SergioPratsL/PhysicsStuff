
Sigma_z = MatrizSpin_4_4(3);

[at, ax, ay, az] = MatricesAlfa();

gamma_t = MatrizGamma(0)

conmut_sigma_z_alfa_x = Sigma_z' * ax - ax' * Sigma_z;

conmut_sigma_z_alfa_y = Sigma_z' * ay - ay' * Sigma_z;

conmut_sigma_z_alfa_z = Sigma_z' * az - az' * Sigma_z;

% V no perjudicará al la conmutación del momento angular.
conmut_sigma_V = Sigma_z' * gamma_t - gamma_t' * Sigma_z;



spinor = [1, 0];
v = [0, 0.4, 0];
p = fGamma(v) * v;

bispinor = DiracSpinorPlainWave(p, spinor)

prod_matrices = ax*p(1) + ay*p(2) + az*p(3);

conmut_H_L = at * prod_matrices;

variacion_onda_Lz = conmut_H_L * bispinor.';

variacion_onda_Lz_horiz = variacion_onda_Lz.'

ratio1 = variacion_onda_Lz(1) / bispinor(1);

ratio4 = variacion_onda_Lz(4) / bispinor(4);

% El componente 4 cambia mucho más rápido.

dif_Sz = conj(bispinor) * Sigma_z * variacion_onda_Lz;


% supongamos la onda plana anterior, [H,Sz] es esto:
conmut_H_Sz = at * (ay*p(1) - ax*p(2));

% [H, Sz]|Phi> vale esto:;
conmut_H_Sz_sobre_Phi = (conmut_H_Sz * bispinor.').'

bisinor_bien_girado = bispinor.';

sz_bispinor = (Sigma_z * bisinor_bien_girado).'


