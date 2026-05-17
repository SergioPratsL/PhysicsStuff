clear;

% El objetivo es, tras ver que la solución de Landau tiene que se buena, y
% tras darme cuenta que con la componnente '4' imaginaria, hay corriente en
% Y, ver que esa corriente sea la adecuada pada el T31-T13.

% Hago algo de trabajo previo para ver en qué direcciones 
% puede haber momento.

% El campo magnético es en Z y Ay=B*x
% La derivada de la onda relevante debería ser sólo en X, incluso aunque
% haya derivada de fase en Y o en Z (creo y espero).

q = 1.602176 * 10^-19;
m_elec = 9.10938 * 10^-31;
perme = 8.854187 * 10^-12;
c = 299792458;
%h = 6.626070 * 10^-34;
%h_bar = h / (2*pi);
h_bar = 6.626070 * 10^-34 / (2*pi);
cte_fina = 1/(4*pi*perme) * q^2 /(h_bar*c);
%permeabilidad = 2 * cte_fina / q^2 * (h/c);
radio_bohr = 4*pi*perme*h_bar^2 / (m_elec*q^2);

mc2 = m_elec * c^2;

[gt, gx, gy, gz] = MatricesGamma();

[at, ax, ay, az] = MatricesAlfa();

[sx, sy, sz] = MatricesPauli();


% Esto representa un spin +Z con parte imaginaria bastante fuerte
onda_base = [1, 0, 0, 1i*0.4].';

% Hacer tres pruebas: 

% Prueba 1, la diferencia de onda es similar a la onda (cambia el módulo).
% Las diferencias en X son imaginarias, las de Y son reales y las de Z son
% nulas.
%dif_onda = -onda_base;

% Prueba 2, la diferencia es en la onda pequeña: cambia k(x).
% Igual: las diferencias en X son imaginarias, las de Y son reales y las de Z son
% nulas.
% dif_onda = -[0,0,0, i].';

% Prueba 3, la diferencia de onda es similar a la onda pero multiplicada
% por i (cambia la fase).
% ¡No cambia nada de nada! Todo es cero, supongo que por los dos términos
% consecuecia de aplicar la derivada en un lado y en otro.
dif_onda = -1i * onda_base;


matriz_dif_px = MatrizFlujoSpin(1, 0, 1);
dif_jx_px = dif_onda' * matriz_dif_px * onda_base + onda_base' * matriz_dif_px * dif_onda;

matriz_dif_py = MatrizFlujoSpin(2, 0, 1);
dif_jy_py = dif_onda' * matriz_dif_py * onda_base + onda_base' * matriz_dif_py * dif_onda;

matriz_dif_pz = MatrizFlujoSpin(3, 0, 1);
dif_jz_pz = dif_onda' * matriz_dif_pz * onda_base + onda_base' * matriz_dif_pz * dif_onda;


% Ahora la parte de probar la función (la de nivel más básico).
B = 1;

x = 2*10^3 * radio_bohr;     % jy_norm = -1.2418e-04    relacion OK,
%x = 10^3 * radio_bohr;       % jy_norm =  -6.2092e-05   relacion OK.
%x = -10^3 * radio_bohr;      % jy_norm =  6.2092e-05   relacion OK.
%x = 5*10^2 * radio_bohr;     % jy_norm = -3.1046e-05   relacion OK.
%x = 10^2 * radio_bohr;        % jy_norm = -6.2092e-06   relacion OK.
%x = 10 * radio_bohr;         % jy_norm = -6.2092e-07   relacion OK.
%x = radio_bohr;              % jy_norm = -6.2092e-08    relacion OK.
% 10^4 es demasiado lejos, gaus es 10^-93!
%x = 10^4 * radio_bohr;        % jy_norm = -6.2092e-04  rel OK

qB = q*B;
hc = h_bar * c;
% NOOOO
% qB_div_h = qB/h;
qB_div_h = qB/h_bar;

E = sqrt(mc2^2 + 2*h_bar*qB*c^2);
factor_peq = -2i*qB*c/(E+mc2);
landau_base = h_bar * qB/m_elec;
landau_dirac = E - mc2;

% Que le den a la gausiana, normalizaré...
gaus = exp(-qB_div_h/2 * x^2);
phi = gaus * [1, 0, 0, factor_peq * x];
dphi_dx = gaus * [-qB_div_h*x, 0, 0, factor_peq*(1 - qB_div_h * x^2)];

% Con estos valores "normalizados", jt deja de tener sentido
%phi = [1, 0, 0, factor_peq * x];
%dphi_dx = [-qB_div_h*x, 0, 0, factor_peq*(1 - qB_div_h * x^2)];

[jt, jx, jy, jz]  = ObtenCorrientesBispinor(phi);

jy_norm = jy / jt

% Absolutamente equivalentes.
py_asociado_a_jy = jy_norm * c * m_elec
dp_y = ObtenDiferenciaFlujoMomentoSoloDx(phi.', dphi_dx.', h_bar) / jt;

% Ignorar el potencial Ay fue la causa que tuviera tantos problemas!
dp_y_corregido = dp_y - qB * x

%py_asociado_a_jy = jy * E / c
%dp_y = ObtenDiferenciaFlujoMomentoSoloDx(phi.', dphi_dx.', h_bar)

ratio = dp_y_corregido / py_asociado_a_jy;

% H1 es dPhi1/dt, por tanto debe ser E*Phi(1)
H1 = 1i*c*h_bar*dphi_dx(4) + 1i*c*qB*x*phi(4) + mc2*phi(1);
H4 = 1i*c*h_bar*dphi_dx(1) - 1i*c*qB*x*phi(1) - mc2*phi(4);

% El Hamiltoniano parece estar bien
dPhi_dt = [H1, 0, 0, H4];
E_H1 = H1 / phi(1) - mc2;
E_H4 = H4 / phi(4) - mc2;
valor_esperable = landau_dirac;


% El completo está en Probaturas_potencial_central_justificar_corrientes.m
function dp = ObtenDiferenciaFlujoMomentoSoloDx(Phi, dPhi_dx, h_bar)
    matriz = MatrizFlujoSpin(2, 0, 1);
    dSuvr_dx = dPhi_dx' * matriz * Phi + Phi' * matriz * dPhi_dx;
    
    % Parece que debe ser multiplicado por h_bar, no h_bar/2 como en el
    % átomo de Hidrógeno (lo cual es una diferencia insufrible).
    dp = dSuvr_dx * h_bar / 2; 
end






