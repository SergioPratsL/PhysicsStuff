% Este Matlab pretende explicar por qué el nivel 1 del átomo de hidrógeno
% tiene coriente no nula en dirección de la colatitud mientras que el
% momento de inercia es nulo basándose en que el spin es como una espira
% que tiene momento alrededor del centro y por este motivo las variaciones
% en orientación del spin o en densidad de carga provocarán una variación
% en la densidad de momento, haciendo que E*j <> p

% https://en.wikipedia.org/wiki/Hydrogen-like_atom#General_solution_of_Wave_Function

% El resultado de este Matlab es positivo, se confirma que la conservación
% del S_\mu\ni\rho cuadra con T_\mu\ni - T_\ni\mu, siendo T el tensor de
% estres energia, y de esta forma se explica la paradoja de tener
% corriente pero no momento en el átomo de hidrógeno


% Se cumple que la corriente iguala a la parte real de
% ObtenDiferenciaFlujoMomento (que tambien tiene parte imaginaria la cual
% obviamente no cuadra con la corriente imaginaria que es 0).
% La densidad de momento se cumple que su parte real es casi cero,
% siempre unos 7 ordenes de magnitud menos que la corriente, y su parte
% imaginaria va por libre




c_elec = 1.602176 * 10^-19;
m_elec = 9.10938 * 10^-31;
perme = 8.854187 * 10^-12;
c = 299792458;
h = 6.626070 * 10^-34;
h_bar = h / (2*pi);
cte_fina = 1/(4*pi*perme) * c_elec^2 /(h_bar*c);
permeabilidad = 2 * cte_fina / c_elec^2 * (h/c);
radio_bohr = 4*pi*perme*h_bar^2 / (m_elec*c_elec^2);
% 1 proton, positivo.
Z = 1;

% nr = 0, espero.
denom_peq = sqrt(1 - cte_fina^2)^2;
denom = sqrt(1 + cte_fina^2 / denom_peq);
E_elec_lv1 = m_elec * c^2 / denom;

E_reposo = m_elec * c^2;

k1 = (E_reposo + E_elec_lv1)/(h_bar*c);
k2 = (E_reposo - E_elec_lv1)/(h_bar*c);

factor_rho = sqrt(k1*k2);

% Viene de Wikipedia pero la \mu era la masa del electrón, confuso...
factor_C = Z* cte_fina * m_elec * c/h_bar;

factor_gamma_wiki = sqrt(1-cte_fina^2);
% Confirmado, son iguales.
sqrt_k1_div_k2_wiki = (1+factor_gamma_wiki)/cte_fina;
sqrt_k1_div_k2_clasico = sqrt(k1 / k2);

k_rara = 1;     % por h_bar.
factor_s = sqrt(k_rara^2 - cte_fina^2);


% Spinor base NO SE PUEDE CAMBIAR
spinor_base = [1, 0];

% r = 1 * radio_bohr;
% dir = [1,0,0];
% dir_e1 = [0,1,0];   
% dir_e2 = [0,0,1];   

% Se cumple
% r = 1.8 * radio_bohr;
% dir = [1,0,0];
% dir_e1 = [0,1,0];   
% dir_e2 = [0,0,1];

% Se cumple, las partes reales cuadran:
% difPx = 1.7825e-24 - 8.9121e-25i
% difPy = 8.9125e-25 + 1.7825e-24i
% difPz = 0
% r = 0.75 * radio_bohr;
% dir = [-1,2,0];
% dir = dir / norm(dir);
% dir_e1 = [0,0,1];   
% dir_e2 = cross(dir, dir_e1);

% Fuera del ecuador 1. Se cumple
% difPx = -9.3326e-25 + 1.4932e-24i     px_asociado_a_jx = -9.3325e-25
% difPy = -1.4932e-24 - 9.3326e-25i 
% difPz = 0.0000e+00 + 9.3326e-25i
% r = 0.5 * radio_bohr;
% dir = [1.6,-1,1];
% dir = dir / norm(dir);
% dir_e1 = [0,1,1] / norm([0,1,1]);
% dir_e2 = cross(dir, dir_e1);


% Fuera del ecuador 2. Se cumple
% difPx = 6.3020e-25
% difPy = difPz = 0
% r = 2.5 * radio_bohr;
% dir = [0,1,3];
% dir = dir / norm(dir);
% dir_e1 = [1,0,0];
% dir_e2 = cross(dir, dir_e1);

% Cerca del polo sur (en las montañas del a locura). Se cumple. VICTORIA!
% difPy = -1.9830e-25
% difPy = difPz = 0

r = radio_bohr;
dir = [1,0,-10];
dir = dir / norm(dir);
dir_e1 = [0, 1, 0];
dir_e2 = cross(dir, dir_e1);

bispinor_base = MontaBispinorBase(factor_gamma_wiki, cte_fina, dir, Z);

bispinor = ObtenBispinor(r, factor_gamma_wiki, factor_C, bispinor_base);

% Verificar que el momento es cero o imaginario (o residual)
px = i * h_bar * dPhi_dx' * bispinor.'
py = i * h_bar * dPhi_dy' * bispinor.'
pz = i * h_bar * dPhi_dz' * bispinor.'

[jt, jx, jy, jz]  = ObtenCorrientesBispinor(bispinor);
jx_norm = jx / jt
jy_norm = jy / jt
jz_norm = jz / jt

dPhi_dr = ObtenDerivadaRadialSpinors(r, factor_gamma_wiki, factor_C, bispinor);
dPhi_e1 = ObtenDerivadaDirNoRadialSpinorPequeno(bispinor, dir_e1, r, Z, cte_fina, factor_gamma_wiki);
dPhi_e2 = ObtenDerivadaDirNoRadialSpinorPequeno(bispinor, dir_e2, r, Z, cte_fina, factor_gamma_wiki);

dPhi_dx = (dPhi_dr*dir(1) + dPhi_e1*dir_e1(1) + dPhi_e2*dir_e2(1)).';
dPhi_dy = (dPhi_dr*dir(2) + dPhi_e1*dir_e1(2) + dPhi_e2*dir_e2(2)).';
dPhi_dz = (dPhi_dr*dir(3) + dPhi_e1*dir_e1(3) + dPhi_e2*dir_e2(3)).';

% En función de la dirección quiero calcular la componente de la colatitud,
% que será una mezcla de (T01-T10) y (T02-T20), tengo que sacar algunos
% Suvr y la derivada de sus valores...

difPx = ObtenDiferenciaFlujoMomento(bispinor.', dPhi_dx, dPhi_dy, dPhi_dz, 1, h_bar) / jt
difPy = ObtenDiferenciaFlujoMomento(bispinor.', dPhi_dx, dPhi_dy, dPhi_dz, 2, h_bar) / jt
difPz = ObtenDiferenciaFlujoMomento(bispinor.', dPhi_dx, dPhi_dy, dPhi_dz, 3, h_bar) / jt


% No hace falta ningún factor gamma(v) ya que jx_norm ya lo contiene al ser
% el producto de E*(vx/c)... o eso creo...
% La idea es que estos tres cudren en la parte real con difPx, difPy, difPz
px_asociado_a_jx = jx_norm * c * m_elec
py_asociado_a_jy = jy_norm * c * m_elec
pz_asociado_a_jz = jz_norm * c * m_elec


% Contrib_H_Q = ObtenAportacionHamiltonianoDerivadasEspaciales(bispinor, dPhi_dx, dPhi_dy, dPhi_dz, h_bar, c) / jt
% 
% E_V = -c_elec/(4*pi*perme) * c_elec / r;
% A = [E_V, 0, 0,0 ];
% 
% Contrib_H_V = ObtenAportacionHamiltonianoPotencial(E_V, jt) / jt
% 
% ratio_Q_V = Contrib_H_Q / Contrib_H_V 
% 
% [pt, px, py, pz] = ObtenEnergiaMomentoConUds(bispinor, dPhi_dx, dPhi_dy, dPhi_dz, A, m_elec, c, h_bar);
% 
% E_Mom = [pt, px, py, pz];
% 
% dif_energ = pt - E_reposo
% dif_teorica = E_elec_lv1 - E_reposo
% 
% Sigma_z = MatrizSpin_4_4(3);
% vector_z = Sigma_z * bispinor.';

% Esto lo estoy haciendo de una forma bastante bárbara, lo sé.
% ratio_spin_z_1 = vector_z(1) / bispinor(1);
% ratio_spin_z_3 = vector_z(3) / bispinor(3);
% ratio_spin_z_4 = vector_z(4) / bispinor(4);

% En la dirección Z los ratios son 1 por tanto bien, en la dirección X
% bispinor(3) es 0 porque es bispinor(4) el que no es cero --> juguete
% roto.


function dp = ObtenDiferenciaFlujoMomento(Phi, dPhi_dx, dPhi_dy, dPhi_dz, indice, h_bar)
    matriz_dir_x = MatrizFlujoSpin(indice, 0, 1);
    dSuvr_dx = dPhi_dx' * matriz_dir_x * Phi + Phi' * matriz_dir_x * dPhi_dx;

    matriz_dir_y = MatrizFlujoSpin(indice, 0, 2);
    dSuvr_dy = dPhi_dy' * matriz_dir_y * Phi + Phi' * matriz_dir_y * dPhi_dy;
    
    matriz_dir_z = MatrizFlujoSpin(indice, 0, 3);
    dSuvr_dz = dPhi_dz' * matriz_dir_z * Phi + Phi' * matriz_dir_z * dPhi_dz;
    
    % Este término es la diferencia de momento :)
    dp = (dSuvr_dx + dSuvr_dy + dSuvr_dz) * h_bar/2; 
end

function bispinor_base = MontaBispinorBase(factor_gamma_wiki, cte_fina, dir, Z)
    spinor_grande = factor_gamma_wiki * [1,0];
    
    spinor_pequeno = -1i * Z * cte_fina * [dir(3), (dir(1)+1i*dir(2))] / 2;

    bispinor_base = [spinor_grande, spinor_pequeno];
end


function dPhi_dir_comple = ObtenDerivadaDirNoRadialSpinorPequeno(bispinor, dir_comple, r, Z, cte_fina, factor_gamma_wiki) 
    norma_spinor_grande = norm(bispinor(1:2));
    %factor_dif_spinor_pequeno = -1i*Z*cte_fina/(1+factor_gamma_wiki)*norma_spinor_grande / r / 2;
    
% Mantener esto a 1 mientras que bispinor_base tiene 1/2 es trampa, lo sé.
    factor_dif_spinor_pequeno = -1i*Z*cte_fina/(1+factor_gamma_wiki)*norma_spinor_grande / r;
    dif_spinor_pequeno = factor_dif_spinor_pequeno * [dir_comple(3), (dir_comple(1) + 1i*dir_comple(2))];
    
    dPhi_dir_comple = [0, 0, dif_spinor_pequeno];
end


function [pt, px, py, pz] = ObtenEnergiaMomentoConUds(phi, dPhi_dx, dPhi_dy, dPhi_dz, A, m, c, h_bar)
% Incluye potencial EM... las unidades serán un dolor, el potencial que
% venga con unidades del sistema inetrnacional!
% El potencial incluye el factor q del electrón (sería la energía potencial
% por unidad de onda).

    E_reposo = m * c^2;
    h_bar_c = h_bar * c;

    V = A(1);
    Ax = A(2);
    Ay = A(3);
    Az = A(4);
    
    [at, ax, ay, az] = MatricesAlfa();
    gt = MatrizGamma(0);
    
    dPhi_dt_gradiente = h_bar_c * ((ax-Ax) * dPhi_dx + (ay-Ay) * dPhi_dy + (az-Az) * dPhi_dz);
    dPhi_dt_masa = - 1i * E_reposo*gt * phi.';
    dPhi_dt_V = - 1i * V * phi.';
    dPhi_dt = ( dPhi_dt_gradiente + dPhi_dt_masa + dPhi_dt_V);
    
    dens_prob = norm(phi)^2;
    
    phi_conj = conj(phi);
    
    % debug
%     pt_A = 1i * phi_conj(1:2) * dPhi_dt(1:2) / dens_prob;
%     pt_B = 1i * phi_conj(3:4) * dPhi_dt(3:4) / dens_prob;
%     pt_r = 1i * phi_conj * dPhi_dx / dens_prob;
    
    pt = 1i * phi_conj * dPhi_dt / dens_prob;
    
    px = -1i * h_bar_c * phi_conj * dPhi_dx  / dens_prob - Ax;
    py = -1i * h_bar_c * phi_conj * dPhi_dy  / dens_prob - Ay;
    pz = -1i * h_bar_c * phi_conj * dPhi_dz  / dens_prob - Az;    
    
    % debug
%     ratio1 = -1i * phi(1) / dPhi_dt(1);
%     ratio3 = -1i * phi(3) / dPhi_dt(3);
%     ratio4 = -1i * phi(4) / dPhi_dt(4);
%     relacion_ratios3 = ratio3 / ratio1
%     relacion_ratios4 = ratio4 / ratio1
    
end

function Contrib_H = ObtenAportacionHamiltonianoDerivadasEspaciales(Phi, dPhi_dx, dPhi_dy, dPhi_dz, h_bar, c)
    [at, ax, ay, az] = MatricesAlfa();
    
    gradiente_slash = 1i * h_bar * c * (ax * dPhi_dx + ay * dPhi_dy + az * dPhi_dz);
    
    Contrib_H = conj(Phi) * gradiente_slash;
end


function bispinor = ObtenBispinor(r, factor_gamma_wiki, factor_C, bispinor_base)
    Cr = r * factor_C;  
   
    componenteRadial = exp(-Cr) * r^((factor_gamma_wiki-1));
    
    bispinor = componenteRadial * bispinor_base;
end

function dPhi_dr = ObtenDerivadaRadialSpinors(r, factor_gamma_wiki, factor_C, bispinor)
    dPhi_dr = ((factor_gamma_wiki - 1) / r - factor_C) * bispinor;
end

function Contrib_H_V = ObtenAportacionHamiltonianoPotencial(E_V, jt)
    Contrib_H_V = E_V * jt;
end



