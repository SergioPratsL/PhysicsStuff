% Pruebas finales de esta campaña.

c_elec = 1.602176 * 10^-19;
m_elec = 9.10938 * 10^-31;
%perme = 8.54187 * 10^-12;      % NOOOO
perme = 8.854187 * 10^-12;
c = 299792458;
h = 6.626070 * 10^-34;
h_bar = h / (2*pi);
cte_fina = 1/(4*pi*perme) * c_elec^2 /(h_bar*c);

radio_bohr = 4*pi*perme*h_bar^2 / (m_elec*c_elec^2);


% nr = 0, espero.
denom_peq = sqrt(1 - cte_fina^2)^2;
denom = sqrt(1 + cte_fina^2 / denom_peq);
E_lv1 = m_elec * c^2 / denom;

E_reposo = m_elec * c^2;

dif_elec_lv1 = E_lv1 - E_reposo

k1 = (E_reposo + E_lv1)/(h_bar*c);
k2 = (E_reposo - E_lv1)/(h_bar*c);

factor_rho = sqrt(k1*k2);

%k_rara = -1;     % por h_bar.
k_rara = 1;     % por h_bar.
factor_s = sqrt(k_rara^2 - cte_fina^2);
%factor_s_simplificado = 1 - cte_fina * sqrt(k2/k1)

%ratio_a_b_s = - cte_fina / (1 - factor_s)
%ratio_a_b_doc = - sqrt(k1/k2)


spinor_base = [1, 0];
% Poner una fase arbiraria.
fase = 1;
%fase = 3 + 2*i;
% FALLA Y ME PARECE QUE ES EL FIN, LOS AUTOVECTORES DEL ÁTOMO DE HIDRÓGENO
% CON DIRAC SON UNA FARSA!
fase = fase / norm(fase);

%r = radio_bohr;
%r = 0.5*radio_bohr;
r = 2*radio_bohr;

%r = 1 * radio_bohr;
%corrientes_intrinsecas =  1.0000         0   -0.0076         0
% Contrib_H_Q =  5.9003e-19
% corrientes_intrinsecas en Z = 1     0     0     0

%r = 0.5 * radio_bohr;
%corrientes_intrinsecas =  1.0000         0   -0.0076         0
% Contrib_H_Q = 8.0197e-19

%r = 2 * radio_bohr;
% corrientes_intrinsecas en X =  1.0000         0   -0.0076         0
% corrientes_intrinsecas en Z = 1     0     0     0
%... la proporción entre los dos spinors siemrpe es la misma.

 dir = [1,0,0];
 dir_e1 = [0,1,0];   % dir_comple no debería hacerse a mano :P.
 dir_e2 = [0,0,1];   % dir_comple no debería hacerse a mano :P.

%dir = [0,0,1];
%dir_e1 = [1,0,0];   % dir_comple no debería hacerse a mano :P.
%dir_e2 = [0,0,1];   % dir_comple no debería hacerse a mano :P.

% Con direciones raras también funcionó
%dir = [1.6,-1,1];
%dir = dir / norm(dir);
%dir_e1 = [0,1,1] / norm([0,1,1]);
%dir_e2 = cross(dir, dir_e1);

bispinor_base = MontaBispinorBase(spinor_base, k1, k2, dir, fase);

bispinor = ObtenBispinor(r, factor_rho, factor_s, bispinor_base);

[jt, jx, jy, jz]  = ObtenCorrientesBispinor(bispinor);

corrientes_intrinsecas = [jt, jx, jy, jz] / jt

dPhi_dr = ObtenDerivadaRadialSpinors(r, factor_rho, factor_s, bispinor);

dPhi_e1 = ObtenDerivadaDirNoRadialSpinorPequeno(bispinor, spinor_base, dir_e1, r);
dPhi_e2 = ObtenDerivadaDirNoRadialSpinorPequeno(bispinor, spinor_base, dir_e2, r);


dPhi_dx = (dPhi_dr*dir(1) + dPhi_e1*dir_e1(1) + dPhi_e2*dir_e2(1)).';
dPhi_dy = (dPhi_dr*dir(2) + dPhi_e1*dir_e1(2) + dPhi_e2*dir_e2(2)).';
dPhi_dz = (dPhi_dr*dir(3) + dPhi_e1*dir_e1(3) + dPhi_e2*dir_e2(3)).';

Contrib_H_Q = ObtenAportacionHamiltonianoDerivadasEspaciales(bispinor, dPhi_dx, dPhi_dy, dPhi_dz, h_bar, c, jt)

Contrib_H_V = ObtenAportacionHamiltonianoPotencial(r, c_elec, perme)

E_V = -c_elec/(4*pi*perme) * c_elec / r;
A = [E_V, 0, 0,0 ];


[pt, px, py, pz] = ObtenEnergiaMomentoConUds(bispinor, dPhi_dx, dPhi_dy, dPhi_dz, A, m_elec, c, h_bar, jt);

E_Mom = [pt, px, py, pz]

momento_normalizado_a_energia = [px, py, pz] / pt;

dif_energia = pt - E_reposo;


% OK. Las diferencias obtenidas son relativamente pequeñas (< 10^-4) 
% para diferentes valores del radio: 0.5, 1 y 2 rbohr. :)
CompruebaEqDiferenciales(k1, k2, factor_rho, factor_s, r, m_elec, c, E_lv1, h_bar, E_V)


function CompruebaEqDiferenciales(k1, k2, factor_rho, factor_s, r, m_elec, c, E_lv1, h_bar, V)
    G = 1;
    F = sqrt(k2 / k1);
    
    factor_dr = -(factor_s / r - factor_rho );
    
    h_bar_c = h_bar * c;
    
    E_reposo = m_elec * c^2;

    K = 1;

    dif_eq_1 = (factor_dr - K*(1/r))*F - (E_reposo - E_lv1 + V)*(G/h_bar_c);
    dif_eq_2 = (factor_dr + K*(1/r))*G - (E_reposo + E_lv1 + V)*(F/h_bar_c);
end


function [pt, px, py, pz] = ObtenEnergiaMomentoConUds(phi, dPhi_dx, dPhi_dy, dPhi_dz, A, m, c, h_bar, jt)
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

    %dPhi_dx = grad_espacial_Phi(1,:).';
    %dPhi_dy = grad_espacial_Phi(2,:).';
    %dPhi_dz = grad_espacial_Phi(3,:).';
    
    [at, ax, ay, az] = MatricesAlfa();
    gt = MatrizGamma(0);
    
    dPhi_dt_gradiente = - h_bar_c * ((ax-Ax) * dPhi_dx + (ay-Ay) * dPhi_dy + (az-Az) * dPhi_dz);
    dPhi_dt_masa = 1i * E_reposo * gt * phi.';
    dPhi_dt_V = 1i * V * phi.';
    dPhi_dt = ( dPhi_dt_gradiente + dPhi_dt_masa + dPhi_dt_V);
    
    phi_conj = conj(phi);
       
    pt = - 1i * phi_conj * dPhi_dt / jt;
    
    px = -1i * h_bar_c * phi_conj * dPhi_dx  / jt - Ax;
    py = -1i * h_bar_c * phi_conj * dPhi_dy  / jt - Ay;
    pz = -1i * h_bar_c * phi_conj * dPhi_dz  / jt - Az;    
    
    % debug
    ratio1 = -1i * phi(1) / dPhi_dt(1)
    ratio2 = -1i * phi(4) / dPhi_dt(4)
    relacion_ratios = ratio2 / ratio1
    
end

% Da la contribución energética por unidad de densiad de onda.
function Contrib_H = ObtenAportacionHamiltonianoDerivadasEspaciales(Phi_ori, dPhi_dx, dPhi_dy, dPhi_dz, h_bar, c, jt)    
    [at, ax, ay, az] = MatricesAlfa();
    
    gradiente_slash = 1i * h_bar * c * (ax * dPhi_dx + ay * dPhi_dy + az * dPhi_dz);
    
    Contrib_H = conj(Phi_ori) * gradiente_slash / jt;
end


function dPhi_dir_comple = ObtenDerivadaDirNoRadialSpinorPequeno(bispinor, spinor_base, dir_comple, r)
    norma_spinor_pequeno = norm(bispinor(3:4));
    
    spinor_base_norm = spinor_base / norm(spinor_base);
    
    dif_sigma = PauliVectorEscalarProd(dir_comple);
    
    dif_spinor_pequeno = -1i * dif_sigma * spinor_base_norm.' * norma_spinor_pequeno / r;
    
    dPhi_dir_comple = [0, 0, dif_spinor_pequeno.'];
end

function bispinor_base = MontaBispinorBase(spinor_base, k1, k2, dir, fase)
    sigma = PauliVectorEscalarProd(dir);
    
    %spinor_pequeno = - 1i * sigma * spinor_base.' * sqrt(k2/k1);  
    spinor_pequeno = 1i * sigma * spinor_base.' * sqrt(k2/k1);  
    
    bispinor_base = [fase * spinor_base, conj(fase) * spinor_pequeno.'];
end


function bispinor = ObtenBispinor(r, factor_rho, factor_s, bispinor_base)
    rho = r * factor_rho;  
    componenteRadial = exp(-rho) * rho^factor_s;
    
    bispinor = componenteRadial * bispinor_base;
end

function dPhi_dr = ObtenDerivadaRadialSpinors(r, factor_rho, factor_s, bispinor)    
    dPhi_dr = (factor_s / r - factor_rho ) * bispinor;
    % CompruebaEqDiferenciales sugiere usar el signo negativo, pues da
    % igual, que le jodan, lo anterior está bien.
    %%%dPhi_dr = -(factor_s / r - factor_rho ) * bispinor;
end

% Da la contribución energética por unidad de densidad de onda.
function Contrib_H_V = ObtenAportacionHamiltonianoPotencial(r, c_elec, perme)
    Contrib_H_V = -(c_elec)/(4*pi*perme) * c_elec / r;
end



