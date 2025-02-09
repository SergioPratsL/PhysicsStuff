
eje_rot = [1, 1.5, -0.8];
eje_rot = eje_rot / norm(eje_rot);

B = ObtenBaseMatrix(eje_rot)

Ident = B * B';     % OK

eje_rot = [1, 0 , 0];
angulo = pi/3;
matrix = RotationMatrixGeneral(eje_rot, angulo);

eje_rot = [1, 1.2 , 0];
eje_rot = eje_rot / norm(eje_rot);
angulo = pi/3;
matrix2 = RotationMatrixGeneral(eje_rot, angulo)


