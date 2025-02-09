function MatrizRotacion = MatrizRotacionSpinDirac(eje, radianes)

M = MatrizRotacionSpin(eje, radianes);

MatrizRotacion = [M, zeros(2); zeros(2), M];

end
