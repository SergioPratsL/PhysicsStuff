% addpath('C:\Fisica\Scripts matlab\HacerTuTeoria\8_origenes')

% Simplemente son dos SRIs que están rotados "por algún motivo" :P
clear;

% B está rotado esto respecto a A
eje_rot = [0,0,1];
angulo = pi /6;

% v_AB = [0.4, 0, 0];
% mat_rot_4_4 = MatrizRotacion_4_4(eje_rot, angulo);
% mat_rot = MatrizRotacion(eje_rot, angulo);
% 
% t = [1,0,0,0];
% 
% % Rot porque aun viene como A
% t_B_rot = Boost(t, v_AB);
% 
% t_B = t_B_rot * mat_rot_4_4;
% v_BA = t_B(2:4) / t_B(1)
% 
% v_AB_rot = v_AB * mat_rot;
% t_B2 = Boost(t, v_AB_rot);
% v_BA2 = t_B2(2:4) / t_B2(1)
% 
% v_BA3 = -v_AB * mat_rot
% v_BA4 = -mat_rot' * v_AB'
% 
% 
% x = [1, 0.7, 2.1, 0.5];
% x_B_rot = Boost(x, v_AB);
% x_B = x_B_rot * mat_rot_4_4
% 
% x_rot = x * mat_rot_4_4;
% x_B2 = Boost(x_rot, v_AB_rot)


v_ab = [0.4, 0, 0];
v_ac = [0, 0.4, 0];

v_bc = Vel_Addition_Law(v_ac, v_ab)
v_cb = Vel_Addition_Law(v_ab, v_ac)


RW = RotacionWigner(v_ab, v_ac);
RW_3 = RW(2:4, 2:4);

v_cb_alt = -v_bc * RW_3
v_bc_alt = -v_cb * RW_3'

x = [1, 1.2, 0.7, 0];
x_c = Boost(x, v_ac)
x_b = Boost(x, v_ab);
x_c2_rot = Boost(x_b, v_bc)
x_c2 = x_c2_rot * RW



