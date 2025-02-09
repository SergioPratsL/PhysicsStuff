function spin_rotado = RotacionSpin( spin, eje, radianes)
% Rota un spinor de Pauli un porron de grados (no tiene que ser
% una rotación diferencial)

eje_norm = eje / norm(eje);

expM = MatrizRotacionSpin(eje_norm, radianes);

spin_rotado = expM * spin';