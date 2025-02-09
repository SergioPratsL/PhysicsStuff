function [d1_V, d1_Ax, d1_Ay, d1_Az] = PrimDerivParcEnh(R_vect, v_vect)
% Package: aceleracion dual

% Generalizacion de PrimDerivParc para cualquier direccion
% R_vect y v_vect deben ser vectores verticales (separados por ';')


% Ojo casos con Rz <> 0 y velocidades fuera del plano XY!!


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

% Az_ = subs(Az);

dV_x = diff(V, 'rx');
dV_y = diff(V, 'ry');
dV_z = diff(V, 'rz');

dAx_x = diff(Ax, 'rx');
dAx_y = diff(Ax, 'ry');
dAx_z = diff(Ax, 'rz');

dAy_x = diff(Ay, 'rx');
dAy_y = diff(Ay, 'ry');
dAy_z = diff(Ay, 'rz');

dAz_x = diff(Az, 'rx');
dAz_y = diff(Az, 'ry');
dAz_z = diff(Az, 'rz');


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


d1_V = [dV_x_, dV_y_, dV_z_];

d1_Ax = [dAx_x_, dAx_y_, dAx_z_];
d1_Ay = [dAy_x_, dAy_y_, dAy_z_];
d1_Az = [dAz_x_, dAz_y_, dAz_z_];

% Mejor no aplicar la direccion
% d1_Ax = dAx_x * dif(1) + dAx_y * dif(2) + dAx_z * dif(3);
% d1_Ay = dAy_x * dif(1) + dAy_y * dif(2) + dAy_z * dif(3);
% d1_Az = dAz_x * dif(1) + dAz_y * dif(2) + dAz_z * dif(3);

% d1_A = [d1_Ax; d1_Ay; d1_Az];



%T = CambioBaseEnh(R_vect, v_vect)

%T_inv = T^(-1);

%R_new_base = T * R_vect;

%dif_new_base = T * dif;

%Rx_new_base = R_new_base(1);
%Ry_new_base = R_new_base(2);
%Rz_new_base = R_new_base(3);


% Si llego al test dual de nivel 2 tendra que incluir Z!! 
%DPot_x_new_base = DerivadasPotencialYinn( Rx_new_base, Ry_new_base, v, 'x');
%DPot_y_new_base = DerivadasPotencialYinn( Rx_new_base, Ry_new_base, v, 'y');
%DPot_z_new_base = DerivadasPotencialYinn( Rx_new_base, Ry_new_base, v, 'z');


%d1_V_new_base = [DPot_x_new_base(1), DPot_y_new_base(1), DPot_z_new_base(1)];

%d1_Ax_new_base = [DPot_x_new_base(2), DPot_y_new_base(2), DPot_z_new_base(2)];

%d1_Ay_new_base = [DPot_x_new_base(3), DPot_y_new_base(3), DPot_z_new_base(3)];

%d1_Az_new_base = [DPot_x_new_base(4), DPot_y_new_base(4), DPot_z_new_base(4)];


%d1_V = dot( d1_V_new_base, dv_new_base);

%d1_Ax_aux = dot( d1_Ax_new_base, dif_new_base);
%d1_Ay_aux = dot( d1_Ax_new_base, dif_new_base);
%d1_Az_aux = dot( d1_Ax_new_base, dif_new_base);


%d1_A_aux = [d1_Ax_aux, d1_Ay_aux, d1_Az_aux];


%d1_A = T_inv * d1_A_aux;




