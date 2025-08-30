function H = HamiltonianoDirac_Prt(Phi, A_tot, dPhi_dx, dPhi_dy, dPhi_dz)
% Respecto al HamiltonianoDirac lo único que hace es reemplazar A con A_tot

m_elec = 9.10938 * 10^-31;
c = 299792458;
h = 6.626070 * 10^-34;
h_bar = h / (2*pi);
q = 1.602176 * 10^-19;


[~, ax, ay, az] = MatricesAlfa();
gt = MatrizGamma(0);

cof_A = q/(c*h_bar);

DPhi_dx = dPhi_dx.' - 1i * cof_A * A_tot{2}*Phi.';
DPhi_dy = dPhi_dy.' - 1i * cof_A * A_tot{3}*Phi.';
DPhi_dz = dPhi_dz.' - 1i * cof_A * A_tot{4}*Phi.';

termino_cinetico = - 1i * h_bar * c * (ax*DPhi_dx + ay*DPhi_dy + az*DPhi_dz);
termino_masa = gt*m_elec*c^2 * Phi.';
termino_V = q * A_tot{1} * Phi.';

H = termino_cinetico + termino_masa + termino_V;
