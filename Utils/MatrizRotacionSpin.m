function MatrizRotacion = MatrizRotacionSpin(eje, radianes)
% A veces reemplazo eje por la velocidad y los radianes por -1i para hacer
% boosts :P

factor = 1i * eje * (radianes/2);

M = PauliVectorEscalarProd(factor);

MatrizRotacion = ExponencialMatriz(M, 0);

end
