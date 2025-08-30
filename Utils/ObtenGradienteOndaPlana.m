function grad_Phi = ObtenGradienteOndaPlana(Phi, p, A)
% función usada con bispinors de Dirac.
% tanto p como A son 4-vectores
% necesita recibir el momento de la partícula, e.g (m_e * v * gamma(v))
% https://quantummechanics.ucsd.edu/ph130a/130_notes/node488.html

    h = 6.626070 * 10^-34;
    h_bar = h / (2*pi);
    c = 299792458;
    q = 1.602176 * 10^-19;

    dPhi_dt = -1i * Phi * (p(1)+q*A(1)) / h_bar;
    % Estos factores (q/c)*A(2) se deben a que para facilitar los boosts,
    % los potenciales que deberían tener un factor 1/c en la parte
    % espacial, no lo tienen, me resisto a hacer boosts con velocidad de la
    % luz no normalizada...
    dPhi_dx = 1i * Phi * (p(2)+(q/c)*A(2)) / h_bar;
    dPhi_dy = 1i * Phi * (p(3)+(q/c)*A(3)) / h_bar;
    dPhi_dz = 1i * Phi * (p(4)+(q/c)*A(4)) / h_bar;
    
    grad_Phi = [dPhi_dt; dPhi_dx; dPhi_dy; dPhi_dz];
end