function [P, M] = DensidadDipolarDirac(Phi)

    [gt, gx, gy, gz] = MatricesGamma();
    Sx = MatrizSpin_4_4(1);
    Sy = MatrizSpin_4_4(2);
    Sz = MatrizSpin_4_4(3);

    Mx = Phi' * gt * Sx * Phi;
    My = Phi' * gt * Sy * Phi;
    Mz = Phi' * gt * Sz * Phi;

    M = [Mx, My, Mz];

    Px = 1i * Phi' * gx * Phi;
    Py = 1i * Phi' * gy * Phi;
    Pz = 1i * Phi' * gz * Phi;

    P = [Px, Py, Pz];
end