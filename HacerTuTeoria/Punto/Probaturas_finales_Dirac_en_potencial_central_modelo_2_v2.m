% Otro enfoque a las pruebsa finales, se basa en la wikipedia:
% https://en.wikipedia.org/wiki/Hydrogen-like_atom#General_solution_of_Wave_Function


c_elec = 1.602176 * 10^-19;
m_elec = 9.10938 * 10^-31;
%perme = 8.54187 * 10^-12;      % NOOOO
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
E_elec_lv1 = m_elec * c^2 / denom

E_reposo = m_elec * c^2;

k1 = (E_reposo + E_elec_lv1)/(h_bar*c);
k2 = (E_reposo - E_elec_lv1)/(h_bar*c);

factor_rho = sqrt(k1*k2);

% La mu era la masa del electrón, nda que ver con la permeabilidad... 
% ¡¡Hijos de la gran puta!!
factor_C = Z* cte_fina * m_elec * c/h_bar;



factor_gamma_wiki = sqrt(1-cte_fina^2);     % Es mi factor_s
% Confirmado, son iguales.
sqrt_k1_div_k2_wiki = (1+factor_gamma_wiki)/cte_fina;
sqrt_k1_div_k2_clasico = sqrt(k1 / k2);

k_rara = 1;     % por h_bar.
factor_s = sqrt(k_rara^2 - cte_fina^2);
%factor_s_simplificado = 1 - cte_fina * sqrt(k2/k1)

%ratio_a_b_s = - cte_fina / (1 - factor_s)
%ratio_a_b_doc = - sqrt(k1/k2)


% Spinor base NO SE PUEDE CAMBIAR
spinor_base = [1, 0];

%r = 0.9 * radio_bohr;

%r = 1 * radio_bohr;

%r = 0.5 * radio_bohr;

r = 2 *radio_bohr;

dir = [1,0,0];
dir_e1 = [0,1,0];   % dir_comple no debería hacerse a mano :P.
dir_e2 = [0,0,1];   % dir_comple no debería hacerse a mano :P.

%dir = [0,0,1];
%dir_e1 = [1,0,0];   % dir_comple no debería hacerse a mano :P.
%dir_e2 = [0,0,1];   % dir_comple no debería hacerse a mano :P.

% Con direciones raras también funcionó
dir = [1.6,-1,1];
dir = dir / norm(dir);
dir_e1 = [0,1,1] / norm([0,1,1]);
dir_e2 = cross(dir, dir_e1);

bispinor_base = MontaBispinorBase(factor_gamma_wiki, cte_fina, dir, Z);

bispinor = ObtenBispinor(r, factor_gamma_wiki, factor_C, bispinor_base);

[jt, jx, jy, jz]  = ObtenCorrientesBispinor(bispinor);

corrientes_intrinsecas = [jt, jx, jy, jz] / jt

dPhi_dr = ObtenDerivadaRadialSpinors(r, factor_gamma_wiki, factor_C, bispinor);
dPhi_e1 = ObtenDerivadaDirNoRadialSpinorPequeno(bispinor, dir_e1, r, Z, cte_fina, factor_gamma_wiki);
dPhi_e2 = ObtenDerivadaDirNoRadialSpinorPequeno(bispinor, dir_e2, r, Z, cte_fina, factor_gamma_wiki);

dPhi_dx = (dPhi_dr*dir(1) + dPhi_e1*dir_e1(1) + dPhi_e2*dir_e2(1)).';
dPhi_dy = (dPhi_dr*dir(2) + dPhi_e1*dir_e1(2) + dPhi_e2*dir_e2(2)).';
dPhi_dz = (dPhi_dr*dir(3) + dPhi_e1*dir_e1(3) + dPhi_e2*dir_e2(3)).';

Contrib_H_Q = ObtenAportacionHamiltonianoDerivadasEspaciales(bispinor, dPhi_dx, dPhi_dy, dPhi_dz, h_bar, c) / jt

E_V = -c_elec/(4*pi*perme) * c_elec / r;
A = [E_V, 0, 0,0 ];

Contrib_H_V = ObtenAportacionHamiltonianoPotencial(E_V, jt) / jt

ratio_Q_V = Contrib_H_Q / Contrib_H_V 

[pt, px, py, pz] = ObtenEnergiaMomentoConUds(bispinor, dPhi_dx, dPhi_dy, dPhi_dz, A, m_elec, c, h_bar);

E_Mom = [pt, px, py, pz];

dif_energ = pt - E_reposo
dif_teorica = E_elec_lv1 - E_reposo

Sigma_z = MatrizSpin_4_4(3);
vector_z = Sigma_z * bispinor.';

% Esto lo estoy haciendo de una forma bastante bárbara, lo sé.
% ratio_spin_z_1 = vector_z(1) / bispinor(1);
% ratio_spin_z_3 = vector_z(3) / bispinor(3);
% ratio_spin_z_4 = vector_z(4) / bispinor(4);

% En la dirección Z los ratios son 1 por tanto bien, en la dirección X
% bispinor(3) es 0 porque es bispinor(4) el que no es cero --> juguete
% roto.


function bispinor_base = MontaBispinorBase(factor_gamma_wiki, cte_fina, dir, Z)
    spinor_grande = factor_gamma_wiki * [1,0];
    
    spinor_pequeno = -1i * Z * cte_fina * [dir(3), (dir(1)+1i*dir(2))] / 2;
    %spinor_pequeno = -1i * Z * cte_fina * [dir(3), (dir(1)+1i*dir(2))];

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
    % No, esto no se puede negociar
    %dPhi_dt_masa = - 1i * E_reposo * phi.';
    dPhi_dt_masa = - 1i * E_reposo*gt * phi.';
    % No, esto no se puede negociar
    %dPhi_dt_V = - 1i * V*gt * phi.';
    dPhi_dt_V = - 1i * V * phi.';
    dPhi_dt = ( dPhi_dt_gradiente + dPhi_dt_masa + dPhi_dt_V);
    
    dens_prob = norm(phi)^2;
    
    phi_conj = conj(phi);
    
    % debug
    pt_A = 1i * phi_conj(1:2) * dPhi_dt(1:2) / dens_prob;
    pt_B = 1i * phi_conj(3:4) * dPhi_dt(3:4) / dens_prob;
    pt_r = 1i * phi_conj * dPhi_dx / dens_prob;
    
    pt = 1i * phi_conj * dPhi_dt / dens_prob;
    
    px = -1i * h_bar_c * phi_conj * dPhi_dx  / dens_prob - Ax;
    py = -1i * h_bar_c * phi_conj * dPhi_dy  / dens_prob - Ay;
    pz = -1i * h_bar_c * phi_conj * dPhi_dz  / dens_prob - Az;    
    
    % debug
    ratio1 = -1i * phi(1) / dPhi_dt(1);
    ratio3 = -1i * phi(3) / dPhi_dt(3);
    ratio4 = -1i * phi(4) / dPhi_dt(4);
    relacion_ratios3 = ratio3 / ratio1
    relacion_ratios4 = ratio4 / ratio1
    
end

function Contrib_H = ObtenAportacionHamiltonianoDerivadasEspaciales(Phi, dPhi_dx, dPhi_dy, dPhi_dz, h_bar, c)
    [at, ax, ay, az] = MatricesAlfa();
    
    gradiente_slash = 1i * h_bar * c * (ax * dPhi_dx + ay * dPhi_dy + az * dPhi_dz);
    
    Contrib_H = conj(Phi) * gradiente_slash;
end


function bispinor = ObtenBispinor(r, factor_gamma_wiki, factor_C, bispinor_base)
    Cr = r * factor_C;  
    
    % Tras la "sucia trampa" hay que restaurar los valores iniciales.
    componenteRadial = exp(-Cr) * r^((factor_gamma_wiki-1));
    %componenteRadial = exp(-2*Cr) * r^(2*(factor_gamma_wiki-1));
    
    bispinor = componenteRadial * bispinor_base;
end

function dPhi_dr = ObtenDerivadaRadialSpinors(r, factor_gamma_wiki, factor_C, bispinor)
    % Tras la "sucia trampa" hay que restaurar los valores iniciales.
    dPhi_dr = ((factor_gamma_wiki - 1) / r - factor_C) * bispinor;
    %dPhi_dr = 2*((factor_gamma_wiki - 1) / r - factor_C) * bispinor;
end

function Contrib_H_V = ObtenAportacionHamiltonianoPotencial(E_V, jt)
    Contrib_H_V = E_V * jt;
end



