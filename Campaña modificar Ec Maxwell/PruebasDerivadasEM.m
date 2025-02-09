

% Prueba uno, campo de Coulomb , partícula se mueve en X o en Y y luego
% hacemos boost a otros sistemas y comparamos dE·E - dB·B
% Por dE o dB se entiend dE/dz y dB/dz siendo dz la 4 velocidad de la
% partícula (es decir vX o vY en estos ejemplos).
% Por simplicidad en el sistema origen "pillamos" a la partícula sobre el 
% eje X de forma que el campo es 1/r^2 (-X)  
% (pongo la carga por delante de la partícula).

% La intensidad inicial es 1.

% Distancia inicial (arbitraria)
r = 1.5;

Ei = [-1/r^2, 0, 0];
Bi = [0,0,0];

% si x crece (por tanto me acerco a la carga)
dE_dx = [2 / r^3, 0, 0];
dE_dy = [0, 1/r^2, 0];
dE_dz = [0, 0, 1/r^2];

dB_dx = [0,0,0];
dB_dy = [0,0,0];
dB_dz = [0,0,0];


v = 0.3;

%Partícula se mueve en X
dEi_dz = v * dE_dx;
dBi_dz = v * dB_dx;

% Pongo punto y coma para mostrar los nuevos resultados
Invariante_mov_x_ori_x = dot(Ei, dEi_dz) - dot(Bi, dBi_dz);

% Boost en X no tiene sentido, comienzo por boost en Y.
v_boost = [0, 0.25, 0];

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost);

% Pongo punto y coma para mostrar los nuevos resultados
Invariante_mov_x_boost_x1 = dot(Ef, dEf_dz) - dot(Bf, dBf_dz);


% Otra pruebecilla        --> OK
v_boost_2 = [0, 0.25, -0.2];

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost_2);
[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost_2);

% Pongo punto y coma para mostrar los nuevos resultados
Invariante_mov_x_boost_x2 = dot(Ef, dEf_dz) - dot(Bf, dBf_dz);


% Y otra                 --> OK
v_boost_3 = [0.18, -0.22, 0.3];

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost_3);
[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost_3);

% Pongo punto y coma para mostrar los nuevos resultados
Invariante_mov_x_boost_x3 = dot(Ef, dEf_dz) - dot(Bf, dBf_dz);



%Partícula se mueve en Y
dEi_dz = v * dE_dy;
dBi_dz = v * dB_dy;

% Pongo punto y coma para mostrar los nuevos resultados
Invariante_mov_x_ori_y = dot(Ei, dEi_dz) - dot(Bi, dBi_dz);


% Boost en X no tiene sentido, comienzo por boost en Y.
v_boost = [0, 0.26, 0];

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost);

% Pongo punto y coma para mostrar los nuevos resultados
Invariante_mov_x_boost_y1 = dot(Ef, dEf_dz) - dot(Bf, dBf_dz);


% Otra pruebecilla        --> OK
v_boost_2 = [0.20, 0, -0.31];

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost_2);
[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost_2);

% Pongo punto y coma para mostrar los nuevos resultados
Invariante_mov_x_boost_y2 = dot(Ef, dEf_dz) - dot(Bf, dBf_dz);


% Y otra                 --> OK
v_boost_3 = [0.10, -0.32, -0.21];

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost_3);
[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost_3);

% Pongo punto y coma para mostrar los nuevos resultados
Invariante_mov_x_boost_y3 = dot(Ef, dEf_dz) - dot(Bf, dBf_dz);


% Pruebas con ondas planas
k = 10;         % c=1 - w=k
w = k;

% Vectores envolventes, no incluyen el coseno(wt-k·x)
E_envolv = [0, 1, 0];
B_envolv = [0, 0, 1];
k_vect = cross(E_envolv, B_envolv);

fase = 0.34 * pi;

% Campos y derivadas dependen de la fase!!
Ei = E_envolv * cos(fase);
dE_dt = - w * E_envolv * sin(fase);
dE_dx = w * E_envolv * sin(fase);
dE_dy = [0, 0, 0];
dE_dz = [0, 0, 0];

Bi = B_envolv * cos(fase);
dB_dt = - w * B_envolv * sin(fase);
dB_dx = w * B_envolv * sin(fase);
dB_dy = [0, 0, 0];
dB_dz = [0, 0, 0];

% Primer bloque de pruebas, partícula mueve en X
v = [0.3, 0, 0];

dEi_dz = v(1) * dE_dx + v(2) * dE_dy + v(3) * v(3) + dE_dt;
dBi_dz = v(1) * dB_dx + v(2) * dB_dy + v(3) * v(3) + dB_dt;


Invariante_mov_x_ori = dot(Ei, dEi_dz) - dot(Bi, dBi_dz)

% Comienzo por boost en Y.
v_boost = [0, 0.26, 0];

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost);

Invariante_mov_x_boost_x_1 = dot(Ef, dEf_dz) - dot(Bf, dBf_dz)

% Un boost raro y carretera (porque el invariante es siempre cero).
v_boost = [-0.5, 0.26, 0.4];

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost);

Invariante_mov_x_boost_x_2 = dot(Ef, dEf_dz) - dot(Bf, dBf_dz)



% Segundo bloque de pruebas, partícula mueve en Y --> Todos dan cero.
v = [0, -0.3, 0];

dEi_dz = v(1) * dE_dx + v(2) * dE_dy + v(3) * v(3) + dE_dt;
dBi_dz = v(1) * dB_dx + v(2) * dB_dy + v(3) * v(3) + dB_dt;

Invariante_mov_y_ori = dot(Ei, dEi_dz) - dot(Bi, dBi_dz)

% Boost 1.
v_boost = [0.3, 0.26, 0];

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost);

Invariante_mov_x_boost_y_1 = dot(Ef, dEf_dz) - dot(Bf, dBf_dz)

% Un boost raro y carretera (porque el invariante es siempre cero).
v_boost = [-0.4, 0.16, -0.3];

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost);

Invariante_mov_x_boost_y_2 = dot(Ef, dEf_dz) - dot(Bf, dBf_dz)




% Tercer bloque de pruebas, partícula mueve en Z 
v = [0, 0, 0.31];

dEi_dz = v(1) * dE_dx + v(2) * dE_dy + v(3) * v(3) + dE_dt;
dBi_dz = v(1) * dB_dx + v(2) * dB_dy + v(3) * v(3) + dB_dt;

Invariante_mov_z_ori = dot(Ei, dEi_dz) - dot(Bi, dBi_dz)


% Boost 1.
v_boost = [0.3, 0.26, 0];

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost);

Invariante_mov_x_boost_z_1 = dot(Ef, dEf_dz) - dot(Bf, dBf_dz)

% Un boost raro y carretera (porque el invariante es siempre cero).
v_boost = [-0.4, 0.16, -0.3];

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost);

Invariante_mov_x_boost_z_2 = dot(Ef, dEf_dz) - dot(Bf, dBf_dz)



% Ahora combinemos a ambas: Coulom mas planas
r = 1.5;

Ei_cou = [-1/r^2, 0, 0];
Bi_cou = [0,0,0];

% si x crece (por tanto me acerco a la carga)
dE_dt_cou = [0,0,0];
dE_dx_cou = [2 / r^3, 0, 0];
dE_dy_cou = [0, 1/r^2, 0];
dE_dz_cou = [0, 0, 1/r^2];

dB_dt_cou = [0,0,0];
dB_dx_cou = [0,0,0];
dB_dy_cou = [0,0,0];
dB_dz_cou = [0,0,0];


k = 4;         % c=1 - w=k
w = k;

% Vectores envolventes, no incluyen el coseno(wt-k·x)
E_envolv_pla = [0, 1, 0];
B_envolv_pla = [0, 0, 1];
k_vect = cross(E_envolv_pla, B_envolv_pla);

fase = 0.34 * pi;

Ei_pla = E_envolv * cos(fase);
dE_dt_pla = - w * E_envolv * sin(fase);
dE_dx_pla = w * E_envolv * sin(fase);
dE_dy_pla = [0, 0, 0];
dE_dz_pla = [0, 0, 0];

Bi_pla = B_envolv * cos(fase);
dB_dt_pla = - w * B_envolv * sin(fase);
dB_dx_pla = w * B_envolv * sin(fase);
dB_dy_pla = [0, 0, 0];
dB_dz_pla = [0, 0, 0];

% Sacamos los campos y las derivadas totales
Ei = Ei_cou + Ei_pla;
Bi = Bi_cou + Bi_pla;

dE_dt = dE_dt_pla + dE_dt_cou;
dE_dx = dE_dx_pla + dE_dx_cou;
dE_dy = dE_dy_pla + dE_dy_cou;
dE_dz = dE_dz_pla + dE_dz_cou;

dB_dt = dB_dt_pla + dB_dt_cou;
dB_dx = dB_dx_pla + dB_dx_cou;
dB_dy = dB_dy_pla + dB_dy_cou;
dB_dz = dB_dz_pla + dB_dz_cou;


% Primer bloque de pruebas, partícula mueve en X
v = [0.3, 0, 0];

dEi_dz = v(1) * dE_dx + v(2) * dE_dy + v(3) * v(3) + dE_dt;
dBi_dz = v(1) * dB_dx + v(2) * dB_dy + v(3) * v(3) + dB_dt;

Invariante_comb_x_ori = dot(Ei, dEi_dz) - dot(Bi, dBi_dz)

% Un boost raro y carretera 
v_boost = [-0.4, 0.16, -0.3];

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost);

Invariante_comb_x_boost_raro = dot(Ef, dEf_dz) - dot(Bf, dBf_dz)


% Segundo bloque de pruebas, partícula mueve en Y
v = [0, -0.3, 0];

dEi_dz = v(1) * dE_dx + v(2) * dE_dy + v(3) * v(3) + dE_dt;
dBi_dz = v(1) * dB_dx + v(2) * dB_dy + v(3) * v(3) + dB_dt;

Invariante_comb_x_ori = dot(Ei, dEi_dz) - dot(Bi, dBi_dz)

% Un boost raro y carretera 
v_boost = [-0.26, -0.25, 0.31];

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost);

Invariante_comb_y_boost_raro = dot(Ef, dEf_dz) - dot(Bf, dBf_dz)



% Tercer bloque de pruebas, partícula mueve en Z
v = [0, 0, 0.3];

dEi_dz = v(1) * dE_dx + v(2) * dE_dy + v(3) * v(3) + dE_dt;
dBi_dz = v(1) * dB_dx + v(2) * dB_dy + v(3) * v(3) + dB_dt;

Invariante_comb_x_ori = dot(Ei, dEi_dz) - dot(Bi, dBi_dz)

% Un boost raro y carretera 
v_boost = [0.1, 0.28, -0.11];

[Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost);

Invariante_comb_z_boost_raro = dot(Ef, dEf_dz) - dot(Bf, dBf_dz)



% PRUEBA(S) DE SUMAR TODAS LAS COMPONENTES A SACO
% Fallo como era de esperar

% Distancia inicial (arbitraria)
r = 1.5;

Ei = [-1/r^2, 0, 0];
Bi = [0,0,0];

% si x crece (por tanto me acerco a la carga)
dE_dx = [2 / r^3, 0, 0];
dE_dy = [0, 1/r^2, 0];
dE_dz = [0, 0, 1/r^2];

dB_dx = [0,0,0];
dB_dy = [0,0,0];
dB_dz = [0,0,0];

v = 0.3;

%Partícula se mueve en X
dEi_dz = v * dE_dx;
dBi_dz = v * dB_dx;
vel = [v, 0, 0];

vi_norm = vel / norm(v);
vi_norm = vi_norm;

% Pongo punto y coma para mostrar los nuevos resultados
%Invariante_Probat_Ori = sum(dEi_dz) - sum(dBi_dz)
Invariante_Probat_Ori = dot(dEi_dz, vi_norm) - dot(dBi_dz, vi_norm)

% Boost en X no tiene sentido, comienzo por boost en Y.
v_boost = [-0.4, 0.25, 0.2];

[dEf_dz, dBf_dz] = Boost_EM(dEi_dz, dBi_dz, v_boost);

vf = Vel_Addition_Law( vel, v_boost);
vf_norm = vf / norm(vf);
vf_norm = vf_norm';

% Pongo punto y coma para mostrar los nuevos resultados
%Invariante_Probat_Boost = sum(dEf_dz) - sum(dBf_dz)
Invariante_Probat_Boost = dot(dEf_dz, vf_norm) - dot(dBf_dz, vf_norm)
