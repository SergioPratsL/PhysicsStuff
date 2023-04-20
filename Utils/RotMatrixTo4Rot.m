
function matrix4 = RotMatrixTo4Rot(matrix3)
% Crea una matriz de tipo [1,0,0,0] y luego 
% para las otras 3 columnas, 0, rotMatrix(1,:)

matrix4_1 = [1,0,0,0];
matrix4_2 = [0, matrix3(1,:)];
matrix4_3 = [0, matrix3(2,:)];
matrix4_4 = [0, matrix3(3,:)];

matrix4 = [matrix4_1; matrix4_2; matrix4_3; matrix4_4];