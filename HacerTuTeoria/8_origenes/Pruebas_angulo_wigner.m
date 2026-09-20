% v_ca = [0.4, 0, 0];
% v_cb = [0, 0.4, 0];


v_ca = [0.6, 0.3, 0.4];
v_cb = [-0.2, 0.4, -0.18];

RW_c = RotacionWigner(v_cb, v_ca);
[dir_rot, angulo] = ExtraeAnguloDeMatriz(RW_c)
%[dir_rot, angulo] = Angulo_Wigner(v_cb, v_ca)

v_ab = Vel_Addition_Law(v_cb, v_ca)

RW_a = RotacionWigner(v_ab, -v_ca);
% Cambiando el orden se invierte la matriz (gira para el otro sentido)
%RW_a = RotacionWigner(-v_ca, v_ab);
[dir_rot2, angulo2] = ExtraeAnguloDeMatriz(RW_a)
%[dir_rot2, angulo2] = Angulo_Wigner(v_ab, -v_ca)

v_ba = Vel_Addition_Law(v_ca, v_cb)

RW_b = RotacionWigner(-v_ba, -v_ca);
[dir_rot3, angulo3] = ExtraeAnguloDeMatriz(RW_b);
% %[dir_rot3, angulo3] = Angulo_Wigner(v_ab_rot, -v_ca)
% [dir_rot3, angulo3] = Angulo_Wigner(-v_ba, -v_ca)

% Angulo2 es el ganador!
dif_2 = angulo - angulo2
dif_3 = angulo - angulo3

% Tras corregir el Angulo_Wigner (comentando mi versión y descomentando la
% de la IA) el angulo3 coincide con el angulo.

% 19.09 Pero dado que el Angulo_Wigner de la IA también demostró estar mal
% en Pruebas_doble_boost_EM.m y RotacionWigner que se basa en otro
% principio lo petó, pues resulta que el angulo que coincide es el angulo2

% Tiene sentido porque v_ab es lo que ve A (calculado desde los datos de C
% pero calculado para A), -v_ca también es lo que ve A que se asume
% alineado con C (cuando calcula C asume que sus datos están alineados con
% los dos).
% Por otro lado -v_ba no es lo que C calcula que ve A sino lo inverso a lo
% que C calcula que ve B, dado que tiene que haber una rotación entre
% ellos, no está bien.
