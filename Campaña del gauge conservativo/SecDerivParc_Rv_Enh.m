function [ddV, ddAx, ddAy, ddAz] = SecDerivParc_Rv_Enh(R_vect, v_vect)
% Package: aceleracion dual

% Generalizacion de la segunda derivada parcial tanto en velocidad como en
% espacio (una derivada en velocidad y otra en espacio).
% R_vect y v_vect deben ser vectores verticales (separados por ';')

% Siempre voy a hacer que en 'diff' la primera derivada sea de velocida y
% la segunda sea espacial... deja de ser simetrica la matriz obtenida!

% Por tanto es muy importante si se ha de trasponer o no!!

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

dV_rx = diff(V, 'rx');
dV_xx = diff(dV_rx, 'vx');
dV_xy = diff(dV_rx, 'vy');
dV_xz = diff(dV_rx, 'vz');

dV_ry = diff(V, 'ry');
dV_yx = diff(dV_ry, 'vx');
dV_yy = diff(dV_ry, 'vy');
dV_yz = diff(dV_ry, 'vz');

dV_rz = diff(V, 'rz');
dV_zx = diff(dV_rz, 'vx');
dV_zy = diff(dV_rz, 'vy');
dV_zz = diff(dV_rz, 'vz');

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


dAx_rx = diff(Ax, 'rx');
dAx_xx = diff(dAx_rx, 'vx');
dAx_xy = diff(dAx_rx, 'vy');
dAx_xz = diff(dAx_rx, 'vz');

dAx_ry = diff(Ax, 'ry');
dAx_yx = diff(dAx_ry, 'vx');
dAx_yy = diff(dAx_ry, 'vy');
dAx_yz = diff(dAx_ry, 'vz');

dAx_rz = diff(Ax, 'rz');
dAx_zx = diff(dAx_rz, 'vx');
dAx_zy = diff(dAx_rz, 'vy');
dAx_zz = diff(dAx_rz, 'vz');

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


dAy_rx = diff(Ay, 'rx');
dAy_xx = diff(dAy_rx, 'vx');
dAy_xy = diff(dAy_rx, 'vy');
dAy_xz = diff(dAy_rx, 'vz');

dAy_ry = diff(Ay, 'ry');
dAy_yx = diff(dAy_ry, 'vx');
dAy_yy = diff(dAy_ry, 'vy');
dAy_yz = diff(dAy_ry, 'vz');

dAy_rz = diff(Ay, 'rz');
dAy_zx = diff(dAy_rz, 'vx');
dAy_zy = diff(dAy_rz, 'vy');
dAy_zz = diff(dAy_rz, 'vz');


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


dAz_rx = diff(Az, 'rx');
dAz_xx = diff(dAz_rx, 'vx');
dAz_xy = diff(dAz_rx, 'vy');
dAz_xz = diff(dAz_rx, 'vz');

dAz_ry = diff(Az, 'ry');
dAz_yx = diff(dAz_ry, 'vx');
dAz_yy = diff(dAz_ry, 'vy');
dAz_yz = diff(dAz_ry, 'vz');

dAz_rz = diff(Az, 'rz');
dAz_zx = diff(dAz_rz, 'vx');
dAz_zy = diff(dAz_rz, 'vy');
dAz_zz = diff(dAz_rz, 'vz');

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
