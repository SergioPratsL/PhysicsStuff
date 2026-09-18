v_ca = [0.4, 0, 0];
v_cb = [0, 0.4, 0];

[dir_rot, angulo] = Angulo_Wigner(v_cb, v_ca)

v_ab = Vel_Addition_Law(v_cb, v_ca)

% El ángulo no cuadra.
[dir_rot2, angulo2] = Angulo_Wigner(v_ab, -v_ca)


mat_rot = MatrizRotacion(dir_rot, angulo);
v_ab_rot = v_ab * mat_rot'

v_ba = Vel_Addition_Law(v_ca, v_cb)

%[dir_rot3, angulo3] = Angulo_Wigner(v_ab_rot, -v_ca)
[dir_rot3, angulo3] = Angulo_Wigner(-v_ba, -v_ca)

dif_2 = angulo - angulo2
dif_3 = angulo - angulo3

% Tras corregir el Angulo_Wigner (comentando mi versión y descomentando la
% de la IA) el angulo3 coincide con el angulo.