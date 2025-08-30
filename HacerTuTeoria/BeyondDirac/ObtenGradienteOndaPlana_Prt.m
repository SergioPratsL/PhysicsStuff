function grad_Phi = ObtenGradienteOndaPlana_Prt(Phi, p, A_tot)
% Version especial de ObtenGradienteOndaPlana con el nuevo término
% A_tot es un array de matrices, la primera vez que lo uso.

    h = 6.626070 * 10^-34;
    h_bar = h / (2*pi);
    c = 299792458;
    q = 1.602176 * 10^-19;

    id = eye(4);

    dPhi_dt = (-1i * (p(1)*id + q*A_tot{1}) * Phi.' / h_bar).';
    % Estos factores (q/c)*{2}) se deben a que para facilitar los boosts,
    % los potenciales que deberían tener un factor 1/c en la parte
    % espacial, no lo tienen, me resisto a hacer boosts con velocidad de la
    % luz no normalizada...
    dPhi_dx = (1i * (p(2)*id + (q/c)*A_tot{2}) * Phi.' / h_bar).';
    dPhi_dy = (1i * (p(3)*id + (q/c)*A_tot{3}) * Phi.' / h_bar).';
    dPhi_dz = (1i * (p(4)*id +  (q/c)*A_tot{4}) * Phi.' / h_bar).';
    
    grad_Phi = [dPhi_dt; dPhi_dx; dPhi_dy; dPhi_dz];
end