% La ideade este script es comprobar que el nuevo término que intento
% introducir se transformará correctamente cuando pase de un SRI a otro,
% Esto implica transformar el tensor E-E con el nuevo término "Prt"

% 09.08.2025. El primer asalto fue un fracaso, ¡no recordaba que la 
% densidad no transforma como un vector! Espero que con el Lagrangiano
% vaya mejor.

clear;

q = 1.602176 * 10^-19;
m_e = 9.10938 * 10^-31;
perm_E = 8.854187 * 10^-12;
c = 299792458;
h = 6.626070 * 10^-34;
h_bar = h / (2*pi);

radio_bohr = 4*pi*perm_E*h_bar^2 / (m_e * q^2);
V_bohr = -1/(4*pi*perm_E) * q / radio_bohr;
E_bohr = - V_bohr / radio_bohr;

[gt, gx, gy, gz] = MatricesGamma();
g5 = gt * gx * gy * gz;

LCM = LeviCivitaMatrix(4);

factor_base = 1 / (m_e * c^2);

spinor_base = [1, 0];

% v es la velocidad de la onda en el SRI original.
%v = [0, 0, 0];
%v = [0.4, 0, 0];
v = [-0.34, 0.41, 0.09];

%u es la velocidad del SRI al que le haremos el boost
%u = [0, -0.4, 0];
%u = [0.4, 0, 0];
u = [-0.09, 0.34, -0.41];

% derivadas de V, de Ax, de Ay y de Az
%D_A = [0,0,0,0; 0,0,0,0; 0,0,0,0; 0,0,0,0]
%D_A = [0,0,0,0; 0,0,0,0; 0,50000,0,0; 0,0,0,0] * E_bohr;
%D_A = [0,50000,100000,200000; 0,0,0,0; 0,0,0,0; 0,0,0,0] * E_bohr;
D_A = [0,50000,100000,200000; 20390,0,-9000,12395; -24853,-118503,0,3000; 6000,-6000,-15000,0] * E_bohr;


% Lo dejaré siempre como 0
%A = [0, 0, 0, 0];
A = [100, -45, 29, 73] * V_bohr;

p = fGamma(v) * v;
% Debería haber un factor c en la energía, pero para los boosts con c=1 es
% mejor no ponerlo.
p4 = m_e * c * [c*fGamma(v), p];

Phi = DiracSpinorPlainWave(p, spinor_base);

A_prt = NuevoTerminoEM_Dirac_Prt(D_A, LCM);

A_tot = {eye(4)*A(1) + A_prt{1}, eye(4)*A(2) + A_prt{2}, eye(4)*A(3) + A_prt{3}, eye(4)*A(4) + A_prt{4}};

grad_Phi = ObtenGradienteOndaPlana_Prt(Phi, p4, A_tot);
dPhi_dx = grad_Phi(2,:);
dPhi_dy = grad_Phi(3,:);
dPhi_dz = grad_Phi(4,:);

T = factor_base * TensorEstresEnergia_QED_Prt(Phi, A_tot, dPhi_dx, dPhi_dy, dPhi_dz)

v_desde_T = T(1, 2:4) / T(1,1)

% Transformación directa.
T_boost = BoostTensorEE(T, u)

% Ahora transformemos todo y evaluémoslo desde el nuevo SRI, es un boost
% "normalizado", por lo tanto no hay que poner un "c" al tiempo.
A_u = Boost(A, u);

v_u = Vel_Addition_Law(v, u)
p_u = fGamma(v_u) * v_u;
p4_u = m_e * c * [c*fGamma(v_u), p_u];

D_A_u = BoostDerivadas4Vector(D_A, u);

Phi_u = DiracSpinorPlainWave(p_u, spinor_base);

A_prt_u = NuevoTerminoEM_Dirac_Prt(D_A_u, LCM);

A_tot_u = {eye(4)*A_u(1) + A_prt_u{1}, eye(4)*A_u(2) + A_prt_u{2}, eye(4)*A_u(3) + A_prt_u{3}, eye(4)*A_u(4) + A_prt_u{4}};

grad_Phi_u = ObtenGradienteOndaPlana_Prt(Phi_u, p4_u, A_tot_u);
dPhi_u_dx = grad_Phi_u(2,:);
dPhi_u_dy = grad_Phi_u(3,:);
dPhi_u_dz = grad_Phi_u(4,:);


T_u = factor_base * TensorEstresEnergia_QED_Prt(Phi_u, A_tot_u, dPhi_u_dx, dPhi_u_dy, dPhi_u_dz)

v_u_desde_T_u = T_u(1, 2:4) / T_u(1,1)
