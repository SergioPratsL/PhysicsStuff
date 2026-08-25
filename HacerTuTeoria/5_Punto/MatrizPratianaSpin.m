function mps = MatrizPratianaSpin(bispinor)
% Prueba para proyectar el spin sobre una matriz

mP_x = MatrizPauli(1);
mP_y = MatrizPauli(2);
mP_z = MatrizPauli(3);


sp_x = mP_x * bispinor';
sp_y = mP_y * bispinor';
sp_z = mP_z * bispinor';

lin_1 = SpinorToVector(sp_x');
lin_2 = SpinorToVector(sp_y');
lin_3 = SpinorToVector(sp_z');

mps = [lin_1; lin_2; lin_3];

end

