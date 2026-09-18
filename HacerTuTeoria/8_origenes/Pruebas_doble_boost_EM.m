clear;

% Las velocidades las sabe B que es el observador, C es el SRIC y A
% representa una partícula que sufrirá fuerza por el campo que en el SRIC
% está en reposo.
%v_ca = [0.4, 0, 0];
%v_cb = [0, 0.4, 0];

% v_ca = [0.4, 0.1, 0];
% v_cb = [0, 0.4, -0.3];
% E_C = [0, 1, 0];
% B_C = [0, 0.4, 0];

v_ca = [0.6, 0.1, 0];
v_cb = [0.15, 0.52, -0.3];
E_C = [3, 1, 0];
B_C = [0, -0.8, 1];

v_ab = Vel_Addition_Law(v_cb, v_ca);
v_ba = Vel_Addition_Law(v_ca, v_cb);
[dir_rot, angulo] = Angulo_Wigner(v_cb, v_ca);

val1 = cross(E_C, B_C);
F_C = TensorEM(E_C, B_C);  
F_B = BoostTensorEE(F_C, v_cb)


%% ESTO CACA.
% F_A = BoostTensorEE(F_C, v_ca);
% Este con -ang y la veloc con +ang
% F_A_rot = RotacionTensor(F_A, dir_rot, -angulo);
% 
% % El campo está definido en C y viaja hasta B por la vía directa y por la
% % indirecta (a través de A).
% 
% mat_rot = MatrizRotacion(dir_rot, angulo);
% 
% v_ab_rot = v_ab * mat_rot;
% %v_ab_rot = v_ab;
% 
% F_B_rot = BoostTensorEE(F_A_rot, v_ab_rot);
% % Esto da mal.
% %F_B_alt2 = BoostTensorEE(F_A, -v_ba);
% 
% 
% F_B_correcto_espero = RotacionTensor(F_B_rot, dir_rot, angulo)
% %F_B_correcto_espero = F_B_rot
% 
% 
% val_espero = norm(F_B_correcto_espero(1, 2:4))
% val = norm(F_B(1, 2:4))
% valB_espero = norm(F_B_correcto_espero(4, 2:4))
% valB = norm(F_B(4, 2:4))
%% ~ESTO CACA.

MB_ca = Tensor_boosts(v_ca);
MB_ab = Tensor_boosts(v_ab);

% Cortesia de Gemini... pues va mal...
F_B_rot = MB_ab * MB_ca * F_C * MB_ca' * MB_ab';

RW = RotacionWigner(v_ca, v_cb);

% Dio bien!!!
%F_B_correcto_espero = RotacionTensor(F_B_rot, dir_rot, angulo)
F_B_correcto = RW' * F_B_rot * RW

valida_correccion = F_B_correcto - F_B

% Comparar las dos matrices de rotacion para ratificar que el metodo
% Angulo_Wigner es una basura
mat_rot_ang_wig = MatrizRotacion_4_4(dir_rot, angulo)
RW = RW

