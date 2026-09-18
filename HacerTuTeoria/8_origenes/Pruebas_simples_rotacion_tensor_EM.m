% Esto no es más que un ejercicio.


%dir_rot = [0, 0, 1];
%angulo = pi /3;
% E = [1, 0, 0];
% B = [0, 0.8, 0];

dir_rot = [2, 0, 1] / norm(dir_rot);
angulo = pi / 5;
E = [1, 0.4, -0.8];
B = [-0.3, 0.5, 0];

mat_rot = MatrizRotacion(dir_rot, angulo);
mat_rot_4_4 = MatrizRotacion_4_4(dir_rot, angulo);



F = TensorEM(E, B);

F_rot1 = mat_rot_4_4 * F * mat_rot_4_4';

E_rot  = mat_rot * E';
B_rot  = mat_rot * B';

F_rot2 = TensorEM(E_rot, B_rot);


