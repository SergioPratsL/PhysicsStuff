ffunction [dRx_dvdr, dRy_dvdr, dRz_dvdr] = SecDeriv_R(R_vect, v_ini)
% Package: velocidad dual
% Derivamos primero respecto la posicion y luego respecto las velocidad

syms rx;
syms ry;
syms rz;

syms vx;
syms vy;
syms vz;

% R es la velocidad desde el SRI acelerado
Rx = (1/(1-vx^2-vy^2-vz^2)^0.5) * (rx - (rx^2+ry^2+rz^2)^(1/2) * vx);
Ry = (1/(1-vx^2-vy^2-vz^2)^0.5) * (ry - (rx^2+ry^2+rz^2)^(1/2) * vy);
Rz = (1/(1-vx^2-vy^2-vz^2)^0.5) * (rz - (rx^2+ry^2+rz^2)^(1/2) * vz);

rx = R_vect(1);
ry = R_vect(2);
rz = R_vect(3);

vx = v_ini(1);
vy = v_ini(2);
vz = v_ini(3);


dRx_x = diff(Rx, 'rx'); 
dRx_xx = diff(dRx_x, 'vx');
dRx_xy = diff(dRx_x, 'vy');
dRx_xz = diff(dRx_x, 'vz');

dRx_y = diff(Rx, 'ry'); 
dRx_yx = diff(dRx_y, 'vx');
dRx_yy = diff(dRx_y, 'vy');
dRx_yz = diff(dRx_y, 'vz');

dRx_z = diff(Rx, 'rz'); 
dRx_zx = diff(dRx_z, 'vx');
dRx_zy = diff(dRx_z, 'vy');
dRx_zz = diff(dRx_z, 'vz');


dRx_xx_ = subs(dRx_xx);
dRx_xy_ = subs(dRx_xy);
dRx_xz_ = subs(dRx_xz);

dRx_yx_ = subs(dRx_yx);
dRx_yy_ = subs(dRx_yy);
dRx_yz_ = subs(dRx_yz);

dRx_zx_ = subs(dRx_zx);
dRx_zy_ = subs(dRx_zy);
dRx_zz_ = subs(dRx_zz);

dRx_dvdr = [dRx_xx_, dRx_xy_, dRx_xz_; dRx_yx_, dRx_yy_, dRx_yz_; dRx_zx_, dRx_zy_, dRx_zz_];


dRy_x =  diff(Ry, 'rx'); 
dRy_xx = diff(dRy_x, 'vx');
dRy_xy = diff(dRy_x, 'vy');
dRy_xz = diff(dRy_x, 'vz');

dRy_y = diff(Ry, 'ry'); 
dRy_yx = diff(dRy_y, 'vx');
dRy_yy = diff(dRy_y, 'vy');
dRy_yz = diff(dRy_y, 'vz');

dRy_z = diff(Ry, 'rz'); 
dRy_zx = diff(dRy_z, 'vx');
dRy_zy = diff(dRy_z, 'vy');
dRy_zz = diff(dRy_z, 'vz');


dRy_xx_ = subs(dRy_xx);
dRy_xy_ = subs(dRy_xy);
dRy_xz_ = subs(dRy_xz);

dRy_yx_ = subs(dRy_yx);
dRy_yy_ = subs(dRy_yy);
dRy_yz_ = subs(dRy_yz);

dRy_zx_ = subs(dRy_zx);
dRy_zy_ = subs(dRy_zy);
dRy_zz_ = subs(dRy_zz);

dRy_dvdr = [dRy_xx_, dRy_xy_, dRy_xz_; dRy_yx_, dRy_yy_, dRy_yz_; dRy_zx_, dRy_zy_, dRy_zz_];



dRz_x =  diff(Rz, 'rx'); 
dRz_xx = diff(dRz_x, 'vx');
dRz_xy = diff(dRz_x, 'vy');
dRz_xz = diff(dRz_x, 'vz');

dRz_y = diff(Rz, 'ry'); 
dRz_yx = diff(dRz_y, 'vx');
dRz_yy = diff(dRz_y, 'vy');
dRz_yz = diff(dRz_y, 'vz');

dRz_z = diff(Rz, 'rz'); 
dRz_zx = diff(dRz_z, 'vx');
dRz_zy = diff(dRz_z, 'vy');
dRz_zz = diff(dRz_z, 'vz');


dRz_xx_ = subs(dRz_xx);
dRz_xy_ = subs(dRz_xy);
dRz_xz_ = subs(dRz_xz);

dRz_yx_ = subs(dRz_yx);
dRz_yy_ = subs(dRz_yy);
dRz_yz_ = subs(dRz_yz);

dRz_zx_ = subs(dRz_zx);
dRz_zy_ = subs(dRz_zy);
dRz_zz_ = subs(dRz_zz);

dRz_dvdr = [dRz_xx_, dRz_xy_, dRz_xz_; dRz_yx_, dRz_yy_, dRz_yz_; dRz_zx_, dRz_zy_, dRz_zz_];

