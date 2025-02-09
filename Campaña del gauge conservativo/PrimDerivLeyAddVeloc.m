function [dS_dvi, dS_dvr] = PrimDerivLeyAddVeloc(v_vect, v_ini)

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


dSx_ux = diff(sx, 'ux'); 
dSx_uy = diff(sx, 'uy'); 
dSx_uz = diff(sx, 'uz'); 

dSy_ux = diff(sy, 'ux'); 
dSy_uy = diff(sy, 'uy'); 
dSy_uz = diff(sy, 'uz'); 

dSz_ux = diff(sz, 'ux'); 
dSz_uy = diff(sz, 'uy'); 
dSz_uz = diff(sz, 'uz'); 

dSx_ux_ = subs(dSx_ux);
dSx_uy_ = subs(dSx_uy);
dSx_uz_ = subs(dSx_uz);

dSy_ux_ = subs(dSy_ux);
dSy_uy_ = subs(dSy_uy);
dSy_uz_ = subs(dSy_uz);

dSz_ux_ = subs(dSz_ux);
dSz_uy_ = subs(dSz_uy);
dSz_uz_ = subs(dSz_uz);


dS_dvi = [dSx_ux_, dSx_uy_, dSx_uz_; dSy_ux_, dSy_uy_, dSy_uz_; dSz_ux_, dSz_uy_, dSz_uz_];



dSx_vx = diff(sx, 'vx'); 
dSx_vy = diff(sx, 'vy'); 
dSx_vz = diff(sx, 'vz'); 

dSy_vx = diff(sy, 'vx'); 
dSy_vy = diff(sy, 'vy'); 
dSy_vz = diff(sy, 'vz'); 

dSz_vx = diff(sz, 'vx'); 
dSz_vy = diff(sz, 'vy'); 
dSz_vz = diff(sz, 'vz'); 

dSx_vx_ = subs(dSx_vx);
dSx_vy_ = subs(dSx_vy);
dSx_vz_ = subs(dSx_vz);

dSy_vx_ = subs(dSy_vx);
dSy_vy_ = subs(dSy_vy);
dSy_vz_ = subs(dSy_vz);

dSz_vx_ = subs(dSz_vx);
dSz_vy_ = subs(dSz_vy);
dSz_vz_ = subs(dSz_vz);


dS_dvr = [dSx_vx_, dSx_vy_, dSx_vz_; dSy_vx_, dSy_vy_, dSy_vz_; dSz_vx_, dSz_vy_, dSz_vz_];