% Lo primero es ver que el TensorEstresEnergia_QED.m no revienta y lo que
% da parece coherente, lo segundo que sin potenciales el tensor transforma
% bien y lo tercero que transforme bien con potenciales.
clear;

m_e = 9.10938 * 10^-31;
c = 299792458;
perm_E = 8.854187 * 10^-12;
q = 1.602176 * 10^-19;
h = 6.626070 * 10^-34;
h_bar = h / (2*pi);

radio_bohr = 4*pi*perm_E*h_bar^2 / (m_e * q^2);
% El menos es porque le casco en negativo a la carga del protón, en vez de
% la del electrón
V_bohr = -1/(4*pi*perm_E) * q / radio_bohr;

factor_base = 1 / (m_e * c^2);

spinor_base = [1, 0];

% Sin potenciales los dos caminos transformaron igual para todas las
% pruebas hechas.

%v = [0, 0, 0];
%v = [0.4, 0, 0];
%v = [0, 0, 0.4];
v = [-0.3, 0.42, 0.4];

%u = [0, 0, 0];
%u = [-0.4, 0, 0];
u = [0.21, -0.29, 0.31];

%A = [0, 0, 0, 0];
%A = [100, 0, 0, 0] * V_bohr;
A = [100, -45, 29, 73] * V_bohr;
% ¡Victoria!


p = fGamma(v) * v;
% Debería haber un factor c en la energía, pero para los boosts con c=1 es
% mejor no ponerlo.
p4 = m_e * c * [c*fGamma(v), p];

Phi = DiracSpinorPlainWave(p, spinor_base);

grad_Phi = ObtenGradienteOndaPlana(Phi, p4, A);
dPhi_dx = grad_Phi(2,:);
dPhi_dy = grad_Phi(3,:);
dPhi_dz = grad_Phi(4,:);

T = factor_base * TensorEstresEnergia_QED(Phi, A, dPhi_dx, dPhi_dy, dPhi_dz)

% Transformación directa.
T_boost = BoostTensorEE(T, u)


% Ahora transformemos todo y evaluémoslo desde el nuevo SRI
A_u = Boost(A, u);

v_u = Vel_Addition_Law(v, u);

p_u = fGamma(v_u) * v_u;
p4_u = m_e * c * [c*fGamma(v_u), p_u];

Phi_u = DiracSpinorPlainWave(p_u, spinor_base);

grad_Phi_u = ObtenGradienteOndaPlana(Phi_u, p4_u, A_u);
dPhi_u_dx = grad_Phi_u(2,:);
dPhi_u_dy = grad_Phi_u(3,:);
dPhi_u_dz = grad_Phi_u(4,:);

T_u = factor_base * TensorEstresEnergia_QED(Phi_u, A_u, dPhi_u_dx, dPhi_u_dy, dPhi_u_dz)
