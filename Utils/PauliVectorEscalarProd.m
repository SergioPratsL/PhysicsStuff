function PV = PauliVectorEscalarProd(V)
% Esta función multiplica un vector por el vector (de matrices) de Pauli

[SigmaX, SigmaY, SigmaZ] = MatricesPauli();

PV = V(1) * SigmaX + V(2) * SigmaY + V(3) * SigmaZ;