% La idea de este script es encontrar las constantes que permiten que 
% ε_μijk γ^i ∂^j A^k nos lleva a tener un termino de momento magnético de
% valor -mu * B, siendo que mu es el momento magnético.
% Recordar que el momento magnético del electrón es:
% -ge*(e/2m)*L, siendo ge el factor anómalo que es casi 2 y siendo L=hbar/2

% Debido a las unidades, pense que se podría añadir un término q/mc
% pero habrá que ver, tengo que hacer los cálculos con todas las unidades.

clear;

c_elec = 1.602176 * 10^-19;
m_elec = 9.10938 * 10^-31;
perme = 8.854187 * 10^-12;
c = 299792458;
h = 6.626070 * 10^-34;
h_bar = h / (2*pi);

[gt, gx, gy, gz] = MatricesGamma();

g5 = gt * gx * gy * gz;

%spinor_base = [1, 0];               % +Z
spinor_base = [1, 1] / sqrt(2);     % +X
%spinor_base = [1, i] / sqrt(2);    % +Y

dir_spin = SpinorToVector(spinor_base)

v = [0, 0, 0];
p = fGamma(v) * v;
E = sqrt(1+norm(p)^2);

phi = DiracSpinorPlainWave(p, spinor_base).';

dens_prob = phi' * phi;
phi = phi / sqrt(dens_prob);

phi_gt = phi' * gt;

% Sz = MatrizSpin_4_4(3);
% Sx = MatrizSpin_4_4(1);
 Sy = MatrizSpin_4_4(2);
% contrib = phi_gt * Sz * phi


contrib2z = phi_gt * g5 * gz * phi;
contrib2x = phi_gt * g5 * gx * phi;
contrib2y = phi_gt * g5 * gy * phi;
contrib2t = phi_gt * g5 * gt * phi;

contribuciones = [contrib2t, contrib2x, contrib2y, contrib2z].';

% matriz_z = 1i * (g5 * gz);
% matriz_x = 1i * g5 * gx;
% matriz_y = g5 * gy;
% matriz_t = g5 * gt;
% 
%  xxx = - 1i * gt * matriz_y;
% % St_obtenido_por_induccion = - 1i * gt * g5 * gt:
% 
% % Anticonmutan
% Gamma_5_Gamma_t_conmut = gt * g5 + g5 * gt;


% La hora de la verdad, esto cambia los autovectores... supongamos que el
% campo es dxAy:
%termino_H = -1i * gt * g5 * gz * phi
termino_H = gt * g5 * gz * phi
phi = phi

dir_spin = SpinorToVector(termino_H(1:2).')

% Una vez le quitas el "i" las contribuciones hacen que precese sobre el
% plano XY, la dir_spin sigue siendo  en sentido opuesto (para X e Y)
% pero el loop siguiente muestra que el bicho precesa :).



% phi_actual = phi;
% factor_dif = 1/1000;
% %matriz_nuevo_termino = -1i * gt * g5 * gz;
% matriz_nuevo_termino = gt * g5 * gz;
% 
% n = 1;
% while n <= 1000
%     termino_H = matriz_nuevo_termino * phi_actual;
%     
%     phi_actual = phi_actual + termino_H / factor_dif;
%     phi_actual = phi_actual / norm(phi_actual);
%     
%     n = n + 1;
% end
% 
% phi_final = phi_actual
% dir_spin_final = SpinorToVector(phi_final(1:2).')
