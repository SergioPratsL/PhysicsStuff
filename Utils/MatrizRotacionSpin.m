function MatrizRotacion = MatrizRotacionSpin(eje, radianes)
% A veces reemplazo eje por la velocidad y los radianes por -1i para hacer
% boosts :P

factor = 1i * eje * (radianes/2);

pl_x = MatrizPauli(1);
pl_y = MatrizPauli(2);
pl_z = MatrizPauli(3);

M = pl_x*factor(1) + pl_y*factor(2) + pl_z*factor(3);

MatrizRotacion = ExponencialMatriz(M, 0);

end
