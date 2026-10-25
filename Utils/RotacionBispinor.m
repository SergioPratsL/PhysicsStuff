function spin_rotado = RotacionBispinor( spin, eje, radianes)
% Rota un spinor de Pauli un porron de grados (no tiene que ser
% una rotación diferencial)

eje_norm = eje / norm(eje);

expM = MatrizRotacionSpinDirac(eje_norm, radianes);

spin_rotado = expM * spin.';