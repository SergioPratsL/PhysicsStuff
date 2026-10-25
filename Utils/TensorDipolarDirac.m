function MP = TensorDipolarDirac(Phi)
% obtiene matriz de densdades dipolares para un spinor de Dirac
% VA MAL :@

    [gt, gx, gy, gz] = MatricesGamma();
    Sx = MatrizSpin_4_4(1);
    Sy = MatrizSpin_4_4(2);
    Sz = MatrizSpin_4_4(3);

    Mx = Phi' * gt * Sx * Phi;
    My = Phi' * gt * Sy * Phi;
    Mz = Phi' * gt * Sz * Phi;

    M = [Mx, My, Mz];

    Px = -1i * Phi' * gx * Phi;
    Py = -1i * Phi' * gy * Phi;
    Pz = -1i * Phi' * gz * Phi;

    P = [Px, Py, Pz];
    
    % El tensor debe llevar -P y M.
    %MP = TensorEM(P, M);
    MP = TensorEM(-P, M);
end