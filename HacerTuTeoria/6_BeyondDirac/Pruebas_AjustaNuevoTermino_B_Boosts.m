% Este script es la continuación del AjustaNuevoTermino_B. Es otro script,
% no es una variación, ¡Son otras pruebas!

clear;

q = 1.602176 * 10^-19;
me = 9.10938 * 10^-31;
perme = 8.854187 * 10^-12;
c = 299792458;
h = 6.626070 * 10^-34;
h_bar = h / (2*pi);

[gt, gx, gy, gz] = MatricesGamma();

g5 = gt * gx * gy * gz;

% En esta prueba el spin será siempre en Z
spinor_base = [1, 0];            % +Z

%v = [0, 0, 0];
%v = [0.4, 0, 0];
v = [0, 0, 0.4];
p = fGamma(v) * v;
E = sqrt(1+norm(p)^2);

phi = DiracSpinorPlainWave(p, spinor_base).';

dens_prob = phi' * phi;
% Para ver cómo transforman los boost, mejor no considerar la densidad de
% probabilidad.
%phi = phi / sqrt(dens_prob);

phi_gt = phi' * gt;
%phi_gt = phi';  Componente t siempre da 0, ¡DESCARTADO!

% Si el componente t siempre da 0, no parece que sea un vector P-L...
% Hay que tener en cuenta que estos cálculosde aqui ignoran el potencial
% de momento sólo examinan cómo se comportan ciertas ondas planas ante las
% matrices...

val_x = -1i * phi_gt * g5 * gx * phi;
val_y = -1i * phi_gt * g5 * gy * phi;
val_z = -1i * phi_gt * g5 * gz * phi;
val_t = -1i * phi_gt * g5 * gt * phi;

resultado = [val_t, val_x, val_y, val_z]

% Así transforma como un vector P-L, debe ser pues, el vector P-L
% resultado con v=0 es [0     0     0    1]          % v = [0, 0, 0]
% resultado = [0         0         0    1.0000]       % v = [0.4, 0, 0]
% resultado = [0.4364         0         0    1.0911]     % v = [0, 0, 0.4]


% Esto era para ver que lo anterior no transforma como la densidad de
% momento dipolar magnético, que transforma como B en el campo EM
%_boost = Boost([0,0,0,1], v);

%[dens_e, dens_m] = Boost_EM([0,0,0], [0, 0, 1], v)
%[dens_m, dens_e] = Boost_EM([0, 0, 1], [0,0,0], v)


% Val z es el término de combinar Gamma_Z con dx y Ay, va a H directamente.
% Esto sería la energía por cada unidad de campo magnético en dirección Z
% expresada en electro-volts, esto aplicaría si tenemos un campo B en Z,
% por ejemplo porque dxAy != 0
cte_ev = val_z * h_bar * (q/me) / q;


% Ver si este término afecta al momento cuando la onda no está en reposo.
% esto sería momento.
contrib_pz = -1i * phi' * gz * g5 * gt * phi
% con v=0 v=0.4X da -1.
% con v=0.4Z da -1.0911




