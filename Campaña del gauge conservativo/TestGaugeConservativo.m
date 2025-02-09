% Sctipt del gauge conservativo
% ¿Es capaz de generar el potencial a partir de sus derivadas parciales?


% Prueba 1. Sale el doble
E = [5, 0, 2];
B = [0, 0.5, 0];


% Prueba 2. Sale el doble.
E = [3, 1, 2];
B = [-1, 0.5, 0.6];


% Prueba 3. Sale el doble.
E = [3, 11, 23];
B = [1, 0, 1];


% Primero la derivada temporal sobre particula quita:
dt_V = 0;
dt_Ax = - E(1);
dt_Ay = - E(2);
dt_Az = - E(3);


% Direccion X:
dx_V = - E(1);
dx_Ax = 0;
dx_Ay = B(3);
dx_Az = - B(2);


% Direccion Y:
dy_V = - E(2);
dy_Ax = - B(3);
dy_Ay = 0;
dy_Az = B(1);


% Direccion Z:
dz_V = - E(3);
dz_Ax = B(2);
dz_Ay = - B(1);
dz_Az = 0;


E_calc_x = - dx_V - dt_Ax;
E_calc_y = - dy_V - dt_Ay;
E_calc_z = - dz_V - dt_Az;

B_calc_x = dy_Az - dz_Ay;
B_calc_y = dz_Ax - dx_Az;
B_calc_z = dx_Ay - dy_Ax;


E = E

E_calc = [E_calc_x, E_calc_y, E_calc_z]

B = B

B_calc = [B_calc_x, B_calc_y, B_calc_z]


