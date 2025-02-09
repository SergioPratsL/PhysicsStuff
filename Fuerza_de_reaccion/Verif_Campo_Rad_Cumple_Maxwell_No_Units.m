% Matlab para verificar que el campo radiado de una carga puntual cumple
% las ecuaciones de Maxwell el solito, sin necesidad de implicar al campo
% inducido (que por tanto las tendrá que cumplir también).

% Quito las unidades porque sospecho que Matlab lo está jodiendo todo!

clear;      

% Prueba 1 - ok
a = (1/1000) * [0, 1, 0];
dt_a = [0,0,0];

% Prueba 2 - ok, el balance no ha de ser cero porque hay dt(a)
a = (1/1000) * [0, 1, 0];
dt_a = 100 * [0,0,1];               


% Prueba 3 - Ok dt, KO los gradientes!!
a = (1/1000) * [1, 0, 0];
dt_a = [0,0,0];               


% Prueba 4 - Ok dt, KO los gradientes!!
%a = (1/1800) * [1, 1, -0.7];
%dt_a = [0,0,0];               

c = 1;

R = [1, 0, 0];      % Distancia 1 metro.

dif = 10^-6;      
dif_t = dif /c;

v = [0, 0, 0];


[E, B] = CampoRadiadoCargaAcelerada_No_Units(R, v, a)


R_dx = R + dif * [1, 0, 0];
[E_dx, B_dx] = CampoRadiadoCargaAcelerada_No_Units(R_dx, v, a);

R_dy = R + dif * [0, 1, 0];
[E_dy, B_dy] = CampoRadiadoCargaAcelerada_No_Units(R_dy, v, a);

R_dz = R + dif * [0, 0, 1];
[E_dz, B_dz] = CampoRadiadoCargaAcelerada_No_Units(R_dz, v, a);

% v_dt será despreciable pero la variación de posición aún lo sería más...
% También asumo que dt_a * dt << a
v_dt = a * dif_t;
a_dt = a + dt_a * dif_t;
[E_dt, B_dt] = CampoRadiadoCargaAcelerada_No_Units(R, v_dt, a_dt);


dx_E = (1/dif)*(E_dx - E);
dx_B = (1/dif)*(B_dx - B);

dy_E = (1/dif)*(E_dy - E);
dy_B = (1/dif)*(B_dy - B);

dz_E = (1/dif)*(E_dz - E);
dz_B = (1/dif)*(B_dz - B);

dt_E = (1/dif_t)*(E_dt - E);
dt_B = (1/dif_t)*(B_dt - B);


dt_Ex_balance = dt_E(1) - c^2*(dy_B(3) - dz_B(2))
dt_Ey_balance = dt_E(2) - c^2*(dz_B(1) - dx_B(3))
dt_Ez_balance = dt_E(3) - c^2*(dx_B(2) - dy_B(1))


dt_Bx_balance = dt_B(1) - (dy_E(3) - dz_E(2))
dt_By_balance = dt_B(2) - (dz_E(1) - dx_E(3))
dt_Bz_balance = dt_B(3) - (dx_E(2) - dy_E(1))

E_Gradient = dx_E(1) + dy_E(2) + dz_E(3)

B_Gradient = dx_B(1) + dy_B(2) + dz_B(3)






