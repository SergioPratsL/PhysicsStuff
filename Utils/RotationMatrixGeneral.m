function matrix = RotationMatrixGeneral( eje_rot, ang )
% Devuelve la matriz para una rotación sobre un eje cualquiera.

val = norm(eje_rot);

if val ~= 0
    eje_rot = eje_rot / val;
end

B = ObtenBaseMatrix(eje_rot);

matrix_base = [1, 0, 0; 0, cos(ang), -sin(ang); 0, sin(ang), cos(ang)];

%matrix = B * matrix_base * B';
matrix = B' * matrix_base * B;

end