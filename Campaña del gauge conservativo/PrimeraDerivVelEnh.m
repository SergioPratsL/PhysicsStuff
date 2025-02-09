function [dV_dv, dAx_dv, dAy_dv, dAz_dv] = PrimeraDerivVelEnh(R_vect, v_vect)
% Package: aceleracion dual

% Los incrementos de velocidad deben llegar ya reducidos, aqui no se
% aplicaran los sigmas


syms rx;
syms ry;
syms rz;

syms vx;
syms vy;
syms vz;

rx = R_vect(1);
ry = R_vect(2);
rz = R_vect(3);

vx = v_vect(1);
vy = v_vect(2);
vz = v_vect(3);


[V, Ax, Ay, Az] = PotYinnFormal;


dV_x = diff(V, 'vx');
dV_y = diff(V, 'vy');
dV_z = diff(V, 'vz');

dAx_x = diff(Ax, 'vx');
dAx_y = diff(Ax, 'vy');
dAx_z = diff(Ax, 'vz');

dAy_x = diff(Ay, 'vx');
dAy_y = diff(Ay, 'vy');
dAy_z = diff(Ay, 'vz');

dAz_x = diff(Az, 'vx');
dAz_y = diff(Az, 'vy');
dAz_z = diff(Az, 'vz');


dV_x_ = subs(dV_x);
dV_y_ = subs(dV_y);
dV_z_ = subs(dV_z);

dAx_x_ = subs(dAx_x);
dAx_y_ = subs(dAx_y);
dAx_z_ = subs(dAx_z);

dAy_x_ = subs(dAy_x);
dAy_y_ = subs(dAy_y);
dAy_z_ = subs(dAy_z);

dAz_x_ = subs(dAz_x);
dAz_y_ = subs(dAz_y);
dAz_z_ = subs(dAz_z);


dV_dv = [dV_x_, dV_y_, dV_z_];

dAx_dv = [dAx_x_, dAx_y_, dAx_z_];
dAy_dv = [dAy_x_, dAy_y_, dAy_z_];
dAz_dv = [dAz_x_, dAz_y_, dAz_z_];





% T = CambioBaseEnh(R_vect, v_vect)
% 
% T_inv = T^(-1);
% 
% R_nb = T * R_vect;
% v_new_base = T * v_vect;
% 
% R = norm(R_vect);
% v = norm(v_vect);
% Sigma = 1 / (1-v^2)^(1/2);
% 
% % Nos olvidamos
% dif_V = 0;
% 
% Denom = (R*v)^2 - dot(R_vect, v_vect)^2
% dAx_dx_1 = (R^2 - R_nb(1)^2)/Denom - 2*v^2 * (R^2 - R_nb(1)^2)^2 / Denom^2;
% 
% Denom = (v * R_nb(2))^2 + (v * R_nb(3))^2;
% Num1 = v^2 * (R - R_nb(1));
% 
% dAx_dx_2 = (-1/v^2) * Num1 / Denom - (1/v) * Num1 * 2 * v * (R_nb(2) + R_nb(3)) / Denom^2;
% 
% dAX_dx_new_base = dAx_dx_1 + dAx_dx_2;
% 
% 
% 
% Denom = (R*v)^2 - dot(R_vect, v_vect)^2
% dAy_dx_1 = (R^2 - R_nb(2)^2)/Denom - 2*v^2*R^2 * (R^2 - R_nb(2)^2)  / Denom^2;
% 
% Denom = (v * R_nb(2))^2 + (v * R_nb(3))^2;
% dAy_dx_2_1 = -1/v^2 * v^2*R_nb(2) / Denom + 1/v * (2*v*R_nb(2)-v*R_nb(1) ) / Denom;
% 
% dAy_dx_2_2 = -1/v * 2 * v * (R_nb(2)^2 + R_nb(3)^2)  * v^2*R_nb(2) / Denom^2;
% 
% 
% dAy_dx_2 = dAy_dx_2_1 + dAy_dx_2_2;
% 
% dAy_dx_new_base = dAy_dx_1 + dAy_dx_2;
% 
