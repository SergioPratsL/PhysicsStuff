function [ddV, ddAx, ddAy, ddAz] = SecDerivParcVelEnh(R_vect, v_vect)
% Package: aceleracion dual

% Generalizacion de SegDerivParc para cualquier direccion de la velocidad
% R_vect y v_vect deben ser vectores verticales (separados por ';')


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

dV_xx = diff(V, 'vx', 2);
dV_x =  diff(V, 'vx');
dV_xy = diff(dV_x, 'vy');
dV_xz = diff(dV_x, 'vz');

dV_yx = dV_xy;
dV_yy = diff(V, 'vy', 2);
dV_y  = diff(V, 'vy');
dV_yz = diff(dV_y, 'vz');

dV_zx = dV_xz;
dV_zy = dV_yz;
dV_zz = diff(V, 'vz', 2);

dV_xx_ = subs(dV_xx);
dV_xy_ = subs(dV_xy);
dV_xz_ = subs(dV_xz);

dV_yx_ = subs(dV_yx);
dV_yy_ = subs(dV_yy);
dV_yz_ = subs(dV_yz);

dV_zx_ = subs(dV_zx);
dV_zy_ = subs(dV_zy);
dV_zz_ = subs(dV_zz);

ddV = [dV_xx_, dV_xy_, dV_xz_; dV_yx_, dV_yy_, dV_yz_; dV_zx_, dV_zy_, dV_zz_];


dAx_xx = diff(Ax, 'vx', 2);
dAx_x  = diff(Ax, 'vx');
dAx_xy = diff(dAx_x, 'vy');
dAx_xz = diff(dAx_x, 'vz');

dAx_yx = dAx_xy;
dAx_yy = diff(Ax, 'vy', 2);
dAx_y =  diff(Ax, 'vy');
dAx_yz = diff(dAx_y, 'vz');

dAx_zx = dAx_xz;
dAx_zy = dAx_yz;
dAx_zz = diff(Ax, 'vz', 2);

dAx_xx_ = subs(dAx_xx);
dAx_xy_ = subs(dAx_xy);
dAx_xz_ = subs(dAx_xz);

dAx_yx_ = subs(dAx_yx);
dAx_yy_ = subs(dAx_yy);
dAx_yz_ = subs(dAx_yz);

dAx_zx_ = subs(dAx_zx);
dAx_zy_ = subs(dAx_zy);
dAx_zz_ = subs(dAx_zz);

ddAx = [dAx_xx_, dAx_xy_, dAx_xz_; dAx_yx_, dAx_yy_, dAx_yz_; dAx_zx_, dAx_zy_, dAx_zz_];


dAy_xx = diff(Ay, 'vx', 2);
dAy_x  = diff(Ay, 'vx');
dAy_xy = diff(dAy_x, 'vy');
dAy_xz = diff(dAy_x, 'vz');

dAy_yx = dAy_xy;
dAy_yy = diff(Ay, vy,2);
dAy_y  = diff(Ay, 'vy');
dAy_yz = diff(dAy_y, 'vz');

dAy_zx = dAy_xz;
dAy_zy = dAy_yz;
dAy_zz = diff(Ay, 'vz', 2);

dAy_xx_ = subs(dAy_xx);
dAy_xy_ = subs(dAy_xy);
dAy_xz_ = subs(dAy_xz);

dAy_yx_ = subs(dAy_yx);
dAy_yy_ = subs(dAy_yy);
dAy_yz_ = subs(dAy_yz);

dAy_zx_ = subs(dAy_zx);
dAy_zy_ = subs(dAy_zy);
dAy_zz_ = subs(dAy_zz);

ddAy = [dAy_xx_, dAy_xy_, dAy_xz_; dAy_yx_, dAy_yy_, dAy_yz_; dAy_zx_, dAy_zy_, dAy_zz_];


dAz_xx = diff(Az, 'vx', 2);
dAz_x  = diff(Az, 'vx');
dAz_xy = diff(dAz_x, 'vy');
dAz_xz = diff(dAz_x, 'vz');

dAz_yx = dAz_xy;
dAz_yy = diff(Az, 'vy',2);
dAz_y  = diff(Az, 'vy');
dAz_yz = diff(dAz_y, 'vz');

dAz_zx = dAz_xz;
dAz_zy = dAz_yz;
dAz_zz = diff(Az, 'vz', 2);

dAz_xx_ = subs(dAz_xx);
dAz_xy_ = subs(dAz_xy);
dAz_xz_ = subs(dAz_xz);

dAz_yx_ = subs(dAz_yx);
dAz_yy_ = subs(dAz_yy);
dAz_yz_ = subs(dAz_yz);

dAz_zx_ = subs(dAz_zx);
dAz_zy_ = subs(dAz_zy);
dAz_zz_ = subs(dAz_zz);

ddAz = [dAz_xx_, dAz_xy_, dAz_xz_; dAz_yx_, dAz_yy_, dAz_yz_; dAz_zx_, dAz_zy_, dAz_zz_];
