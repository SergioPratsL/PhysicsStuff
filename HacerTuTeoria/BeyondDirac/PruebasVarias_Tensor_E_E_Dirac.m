clear;

[gt, gx, gy, gz] = MatricesGamma();

spinor_base = [1, 0];

%v = [0, 0, 0.4];
v = [0.4, 0, 0];
p = fGamma(v) * v;
E = sqrt(1+norm(p)^2);

phi = DiracSpinorPlainWave(p, spinor_base).';

dens_prob = phi' * phi;
phi = phi / sqrt(dens_prob);

%sigma_zz = norm(p) * phi' * gt * gz * phi;
%sigma_zx = norm(p) * phi' * gt * gx * phi;

phi_gt = phi' * gt;


jx_control = phi_gt * gx * phi

E_control = phi' * (gt*p(1)*gx + gt) * phi

E_phi_control = (gt*p(1)*gx + gt) * phi

producto_espero_de_px = phi' * gt * gx * (gt*p(1)*gx + gt) * phi

termino_gt = phi' * gx * phi


% OK, la parte real es 0, la imaginaria no pero eso da igual.
control_1 = phi' * gt *gy * gt * gx * phi;
control_2 = phi' * gt *gy * gt * phi;



v2 = [0.25, 0, 0.32];
p2 = fGamma(v2) * v2;
E2 = sqrt(1+norm(p2)^2);

phi2 = DiracSpinorPlainWave(p2, spinor_base).';

dens_prob2 = phi2' * phi2;
phi2 = phi2 / sqrt(dens_prob2);


% Evaluando en otro punto
pos_x = 4.5;
phi = phi * exp(1i*2*pi*pos_x*v(1));
phi2 = phi2 * exp(1i*2*pi*pos_x*v2(1));

dens_prob_mixta = norm(phi+phi2)^2;
px_mixto = (phi' + phi2') * (v(1)*phi + v2(1)*phi2) / dens_prob_mixta;

jx_mixta = (phi' + phi2') * gt * gx * (phi + phi2) / dens_prob_mixta;
%a = phi2' * gt * gx * phi / dens_prob_mixta
%b = phi' * gt * gx * phi2 / dens_prob_mixta
%c = phi' * gt * gx * phi / dens_prob_mixta
%d = phi2' * gt * gx * phi2 / dens_prob_mixta


coeficiente_V = 
