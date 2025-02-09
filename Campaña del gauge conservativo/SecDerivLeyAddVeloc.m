function [dSx_dv2, dSy_dv2, dSz_dv2] = SecDerivLeyAddVeloc(v_vect, v_ini)
% Package: velocidad dual

% Da las segundas derivadas a la transformacion de velocidades respecto a
% la velocidad del destino (la v) y no la u que es el origen

% Al usar directamente la formula de adicion se elimina el misticismo del
% que tiraba, la matriz "M" y otras mierdas :D.

syms ux;
syms uy;
syms uz;

syms vx;
syms vy;
syms vz;

% 's' es como se denomina a la velocidad compuesta
sx = (vx - ux)/(1 - ux*vx - uy*vy - uz*vz);
sy = (vy - uy)/(1 - ux*vx - uy*vy - uz*vz);
sz = (vz - uz)/(1 - ux*vx - uy*vy - uz*vz);

ux = v_ini(1);
uy = v_ini(2);
uz = v_ini(3);

vx = v_vect(1);
vy = v_vect(2);
vz = v_vect(3);


dSx_x = diff(sx, 'vx'); 
dSx_xx = diff(dSx_x, 'ux');
dSx_xy = diff(dSx_x, 'uy');
dSx_xz = diff(dSx_x, 'uz');

dSx_y = diff(sx, 'vy');
dSx_yx = diff(dSx_y, 'ux');
dSx_yy = diff(dSx_y, 'uy');
dSx_yz = diff(dSx_y, 'uz');

dSx_z = diff(sx, 'vz');
dSx_zx = diff(dSx_z, 'ux');
dSx_zy = diff(dSx_z, 'uy');
dSx_zz = diff(dSx_z, 'uz');


dSx_xx_ = subs(dSx_xx);
dSx_xy_ = subs(dSx_xy);
dSx_xz_ = subs(dSx_xz);

dSx_yx_ = subs(dSx_yx);
dSx_yy_ = subs(dSx_yy);
dSx_yz_ = subs(dSx_yz);

dSx_zx_ = subs(dSx_zx);
dSx_zy_ = subs(dSx_zy);
dSx_zz_ = subs(dSx_zz);

dSx_dv2 = [dSx_xx_, dSx_xy_, dSx_xz_; dSx_yx_, dSx_yy_, dSx_yz_; dSx_zx_, dSx_zy_, dSx_zz_];



dSy_x = diff(sy, 'vx');
dSy_xx = diff(dSy_x, 'ux');
dSy_xy = diff(dSy_x, 'uy');
dSy_xz = diff(dSy_x, 'uz');

dSy_y = diff(sy, 'vy');
dSy_yx = diff(dSy_y, 'ux');
dSy_yy = diff(dSy_y, 'uy');
dSy_yz = diff(dSy_y, 'uz');

dSy_z = diff(sy, 'vz');
dSy_zx = diff(dSy_z, 'ux');
dSy_zy = diff(dSy_z, 'uy');
dSy_zz = diff(dSy_z, 'uz');


dSy_xx_ = subs(dSy_xx);
dSy_xy_ = subs(dSy_xy);
dSy_xz_ = subs(dSy_xz);

dSy_yx_ = subs(dSy_yx);
dSy_yy_ = subs(dSy_yy);
dSy_yz_ = subs(dSy_yz);

dSy_zx_ = subs(dSy_zx);
dSy_zy_ = subs(dSy_zy);
dSy_zz_ = subs(dSy_zz);

dSy_dv2 = [dSy_xx_, dSy_xy_, dSy_xz_; dSy_yx_, dSy_yy_, dSy_yz_; dSy_zx_, dSy_zy_, dSy_zz_];



dSz_x = diff(sz, 'vx');
dSz_xx = diff(dSz_x, 'ux');
dSz_xy = diff(dSz_x, 'uy');
dSz_xz = diff(dSz_x, 'uz');

dSz_y = diff(sz, 'vy');
dSz_yx = diff(dSz_y, 'ux');
dSz_yy = diff(dSz_y, 'uy');
dSz_yz = diff(dSz_y, 'uz');

dSz_z = diff(sz, 'vz');
dSz_zx = diff(dSz_z, 'ux');
dSz_zy = diff(dSz_z, 'uy');
dSz_zz = diff(dSz_z, 'uz');


dSz_xx_ = subs(dSz_xx);
dSz_xy_ = subs(dSz_xy);
dSz_xz_ = subs(dSz_xz);

dSz_yx_ = subs(dSz_yx);
dSz_yy_ = subs(dSz_yy);
dSz_yz_ = subs(dSz_yz);

dSz_zx_ = subs(dSz_zx);
dSz_zy_ = subs(dSz_zy);
dSz_zz_ = subs(dSz_zz);

dSz_dv2 = [dSz_xx_, dSz_xy_, dSz_xz_; dSz_yx_, dSz_yy_, dSz_yz_; dSz_zx_, dSz_zy_, dSz_zz_];


% dSx_x_ = subs(dSx_x)
% dSx_y_ = subs(dSx_y)
% dSx_z_ = subs(dSx_z)
% 
% dSy_x_ = subs(dSy_x)
% dSy_y_ = subs(dSy_y)
% dSy_z_ = subs(dSy_z)
% 
% dSz_x_ = subs(dSz_x)
% dSz_y_ = subs(dSz_y)
% dSz_z_ = subs(dSz_z)

