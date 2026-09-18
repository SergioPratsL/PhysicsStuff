% Estas pruebas son para ver cómo cambian las diferencias de momento
% Se evalúa como cambiará dp al aplicar estas horribles rotaciones.

% R y m están aquí de manera testimonial
R = 1;
m = 1;

v = 0.4;
fase = 0;
gamma = fGamma(v);

% Mi definición de infinitesimal
dt_prim = 10^-4;

% Frecuecia angular en el SRI local de la partícula en MCU
w_prim = gamma * v / R;

% Aceleración que se ve en el SRI local...
a_abs = w_prim^2 * R;

v_CA = v * [-sin(fase), cos(fase), 0];
a = -a_abs * [cos(fase), sin(fase), 0];

p_CA = gamma * [1, v_CA];

% Cambio de 4-momento en el SRI local por un instante dt':
dp_prim = [0, a * dt_prim];

% Este es el valor de toda la vida que en el MCU está bien.
dp_CA_esperado = Boost(dp_prim, v_CA)



[dir_rot, angulo_W_T] = Angulo_Wigner(-v_CA, dp_prim(2:4));
mat_rot4 = MatrizRotacion_4_4(dir_rot, angulo_W_T);
mat_rot_inv = MatrizRotacion(dir_rot, -angulo_W_T);
mat_rot_inv4 = MatrizRotacion_4_4(dir_rot, -angulo_W_T);

dp_prim_rot = dp_prim * mat_rot_inv4;
v_CA_rot = v_CA * mat_rot_inv;

dp_rot = Boost(dp_prim_rot, v_CA_rot);  
% Este es el resultado de dv.
dp = dp_rot * mat_rot4




