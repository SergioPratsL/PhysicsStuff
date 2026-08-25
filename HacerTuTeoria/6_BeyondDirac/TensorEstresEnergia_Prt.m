function T = TensorEstresEnergia_QED(Phi, A, dPhi_dx, dPhi_dy, dPhi_dz)
% Se llama "QED" pero también sirve para Dirac (entiendo que Dirac seria
% sin potenciales y QED con potenciales).
% Este bicho generará el tensor de estres energía para un Dirac  con
% potenciales... lo llaman "Lagrangiano QED"

% https://en.wikipedia.org/wiki/Lagrangian_(field_theory)#Quantum_electrodynamic_Lagrangian
% https://physics.stackexchange.com/questions/414556/energy-momentum-tensor-of-transformed-dirac-lagrangian

c = 299792458;
q = 1.602176 * 10^-19;
h = 6.626070 * 10^-34;
h_bar = h / (2*pi);

T = zeros(4);

[gt, gx, gy, gz] = MatricesGamma();
Phi_gt = conj(Phi * gt);

H = HamiltonianoDirac(Phi, A, dPhi_dx, dPhi_dy, dPhi_dz);

DPhi_dt = H - A(1) * q * Phi.';

DPhi_dx = - 1i * c * (h_bar * dPhi_dx - 1i*(q/c)*A(2)*Phi).';
DPhi_dy = - 1i * c * (h_bar * dPhi_dy - 1i*(q/c)*A(3)*Phi).';
DPhi_dz = - 1i * c * (h_bar * dPhi_dz - 1i*(q/c)*A(4)*Phi).';

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

%dens_prob = Phi * Phi';
%T = T / dens_prob;

end