function T = TensorEstresEnergia_QED_Prt(Phi, A_tot, dPhi_dx, dPhi_dy, dPhi_dz)
% El potencial A ha sido modificado por A_tot, que es un vector de
% matrices.

% https://en.wikipedia.org/wiki/Lagrangian_(field_theory)#Quantum_electrodynamic_Lagrangian
% https://physics.stackexchange.com/questions/414556/energy-momentum-tensor-of-transformed-dirac-lagrangian

c = 299792458;
q = 1.602176 * 10^-19;
h = 6.626070 * 10^-34;
h_bar = h / (2*pi);

T = zeros(4);

[gt, gx, gy, gz] = MatricesGamma();
Phi_gt = conj(Phi * gt);

H = HamiltonianoDirac_Prt(Phi, A_tot, dPhi_dx, dPhi_dy, dPhi_dz);

DPhi_dt = H - q * A_tot{1} * Phi.';

DPhi_dx = - 1i * c * (h_bar * dPhi_dx.' - 1i*(q/c)*A_tot{2}*Phi.');
DPhi_dy = - 1i * c * (h_bar * dPhi_dy.' - 1i*(q/c)*A_tot{3}*Phi.');
DPhi_dz = - 1i * c * (h_bar * dPhi_dz.' - 1i*(q/c)*A_tot{4}*Phi.');

T(1,1) = conj(Phi) * DPhi_dt;

T(1,2) = Phi_gt * gt * DPhi_dx;
T(1,3) = Phi_gt * gt * DPhi_dy;
T(1,4) = Phi_gt * gt * DPhi_dz;

T(2,1) = Phi_gt * gx * DPhi_dt;
T(3,1) = Phi_gt * gy * DPhi_dt;
T(4,1) = Phi_gt * gz * DPhi_dt;

T(2,2) = Phi_gt * gx * DPhi_dx;
T(2,3) = Phi_gt * gx * DPhi_dy;
T(2,4) = Phi_gt * gx * DPhi_dz;

T(3,2) = Phi_gt * gy * DPhi_dx;
T(3,3) = Phi_gt * gy * DPhi_dy;
T(3,4) = Phi_gt * gy * DPhi_dz;

T(4,2) = Phi_gt * gz * DPhi_dx;
T(4,3) = Phi_gt * gz * DPhi_dy;
T(4,4) = Phi_gt * gz * DPhi_dz;

end