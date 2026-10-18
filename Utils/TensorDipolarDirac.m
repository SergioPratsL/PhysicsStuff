function MP = TensorDipolarDirac(Phi)
% obtiene matriz de densdades dipolares para un spinor de Dirac

    [gt, gx, gy, gz] = MatricesGamma();
    Sx = MatrizSpin_4_4(1);
    Sy = MatrizSpin_4_4(2);
    Sz = MatrizSpin_4_4(3);

    Mx = Phi' * Sx * Phi;
    My = Phi' * Sy * Phi;
    Mz = Phi' * Sz * Phi;

    M = [Mx, My, Mz];

    Px = - 1i * Phi' * gx * Phi;
    Py = - 1i * Phi' * gy * Phi;
    Pz = - 1i * Phi' * gz * Phi;

    P = [Px, Py, Pz];

    % Es necesario usar -P para que los tensores tengan la misma formma, según
    % Gemini.
    MP = TensorEM(-P, M);

end