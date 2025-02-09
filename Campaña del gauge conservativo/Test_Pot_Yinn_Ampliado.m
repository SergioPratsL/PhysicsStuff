
% phi indica la direccion de la velocidad en el plano XY

% Prueba 1. OK
Rx = 0;
Ry = 1;
v = 0.4;
phi = 0.00 * pi;


% Prueba 2. OK
Rx = 0;
Ry = 1;
v = 0.4;
phi = 0.29 * pi;


% Prueba 3. OK
Rx = 0.6;
Ry = 0.9;
v = 0.4;
phi = -0.20 * pi;


% Prueba 4. OK
Rx = -1;
Ry = 1.9;
v = 0.428;
phi = 0.32 * pi;


% Prueba 5. OK
Rx = -1;
Ry = -1.9;
v = 0.25;
phi = -0.12 * pi;


% Prueba 6. OK
Rx = 4
Ry = -0.5;
v = -0.14;
phi = 0.17 * pi;


v_vectorial = v * [cos(phi), sin(phi), 0];
R = [Rx, Ry, 0];

% Primero  saco el potencial de Yinn de la forma clasica sobre un vector
% rotado y luego desroto los potenciales

% La rotacion es en el sentido inverso
phi2 = -phi;

Rx_rot = Rx * cos(phi2) - Ry * sin(phi2);
Ry_rot = Ry * cos(phi2) + Rx * sin(phi2);

R_rot = [Rx_rot, Ry_rot, 0];
v_vect = [v, 0, 0];
[V1, A1] = PotencialDeYinn( R_rot, v_vect);

% Desrotar:
A1_desrot(1) = A1(1) * cos(phi) - A1(2) * sin(phi);
A1_desrot(2) = A1(2) * cos(phi) + A1(1) * sin(phi);
A1_desrot(3) = 0;


% Obtener el potencial por la formula directa:
[V2, A2] = PotencialDeYinnAmpliado(R, v_vectorial)


Dif_x = A1_desrot(1) - A2(1)

Dif_y = A1_desrot(2) - A2(2)

Dif_z = A1_desrot(3) - A2(3)