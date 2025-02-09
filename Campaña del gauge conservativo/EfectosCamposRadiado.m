function Dif_Pot_EM = EfectosCamposRadiado(R_vect, v_vect, despl, ac_emis, dv_recept)
% Package: aceleracion dual

% Primera derivada de la variacion de energia y momento del receptor debida
% al campo electromagnetico radiado generado por un instante de aceleracion
% del emisor, al centrarse en "un instante del emisor" y no del receptor
% tiene Doppler.

% Se evalua desde el SRI del emisor
v_ini = [0,0,0];

syms rx;
syms ry;
syms rz;

syms vx;
syms vy;
syms vz;

syms ax; 
syms ay;
syms az;

% Doppler
syms v2x;
syms v2y;
syms v2z;


Dop_inv = (rx^2+ry^2+rz^2)^(1/2) / ((rx^2+ry^2+rz^2)^(1/2) -( rx*v2x + ry*v2y + rz*v2z));

dDop_inv_rx =  diff(Dop_inv, 'rx'); 
dDop_inv_ry =  diff(Dop_inv, 'ry');
dDop_inv_rz =  diff(Dop_inv, 'rz');

dDop_inv_v2x =  diff(Dop_inv, 'v2x'); 
dDop_inv_v2y =  diff(Dop_inv, 'v2y');
dDop_inv_v2z =  diff(Dop_inv, 'v2z');

dDop_inv_rx_ = subs(dDop_inv_rx);
dDop_inv_ry_ = subs(dDop_inv_ry);
dDop_inv_rz_ = subs(dDop_inv_rz);

dDop_inv_v2x_ = subs(dDop_inv_v2x);
dDop_inv_v2y_ = subs(dDop_inv_v2y);
dDop_inv_v2z_ = subs(dDop_inv_v2z);

dDoppler_inv = dDop_inv_rx_ * despl(1) + dDop_inv_ry_ * despl(2) + dDop_inv_rz_ * despl(3);
dDoppler_inv = dDoppler_inv + dDop_inv_v2x_ * dv_recept(1) + dDop_inv_v2y_ * dv_recept(2) + dDop_inv_v2z_ * dv_recept(3);


Ex = (-ax * (rx^2+ry^2+rz^2 - (rx^2+ry^2+rz^2)^(1/2)*(rx*vx+ry*vy+rz*vz)) + ( rx + vx*(rx^2+ry^2+rz^2)^(1/2) ) * (rx*ax+ry*ay+rz*az)) / ( (rx^2+ry^2+rz^2)^(1/2) -(rx*vx+ry*vy+rz*vz) )^3;

Ey = (-ay * (rx^2+ry^2+rz^2 - (rx^2+ry^2+rz^2)^(1/2)*(rx*vx+ry*vy+rz*vz)) + (ry + vy*(rx^2+ry^2+rz^2)^(1/2)) * (rx*ax+ry*ay+rz*az)) / ((rx^2+ry^2+rz^2)^(1/2) -(rx*vx+ry*vy+rz*vz))^3;

Ez = (-az * (rx^2+ry^2+rz^2 - (rx^2+ry^2+rz^2)^(1/2)*(rx*vx+ry*vy+rz*vz)) + (rz + vz*(rx^2+ry^2+rz^2)^(1/2)) * (rx*ax+ry*ay+rz*az)) / ((rx^2+ry^2+rz^2)^(1/2) -(rx*vx+ry*vy+rz*vz))^3;


Bx = ((rz*ay - ry*az) * (rx^2+ry^2+rz^2 - (rx^2+ry^2+rz^2)^(1/2)*(rx*vx+ry*vy+rz*vz)) + (ry*vz - rz*vy) * (rx*ax+ry*ay+rz*az)) / ((rx^2+ry^2+rz^2)^(1/2) -(rx*vx+ry*vy+rz*vz))^3 /(rx^2+ry^2+rz^2)^(1/2);

By = ((rx*az - rz*ax) * (rx^2+ry^2+rz^2 - (rx^2+ry^2+rz^2)^(1/2)*(rx*vx+ry*vy+rz*vz)) + (rz*vx - rx*vz) * (rx*ax+ry*ay+rz*az)) / ((rx^2+ry^2+rz^2)^(1/2) -(rx*vx+ry*vy+rz*vz))^3 /(rx^2+ry^2+rz^2)^(1/2);

Bz = ((ry*ax - rx*ay) * (rx^2+ry^2+rz^2 - (rx^2+ry^2+rz^2)^(1/2)*(rx*vx+ry*vy+rz*vz)) + (rx*vy - ry*vx) * (rx*ax+ry*ay+rz*az)) / ((rx^2+ry^2+rz^2)^(1/2) -(rx*vx+ry*vy+rz*vz))^3 /(rx^2+ry^2+rz^2)^(1/2);


rx = R_vect(1);
ry = R_vect(2);
rz = R_vect(3);

vx = v_ini(1);
vy = v_ini(2);
vz = v_ini(3);

ax = ac_emis(1);
ay = ac_emis(2);
az = ac_emis(3);

v2x = v_vect(1);
v2y = v_vect(2);
v2z = v_vect(3);

dDoppler_inv_ = subs(dDoppler_inv);


dEx_rx =  diff(Ex, 'rx'); 
dEx_ry =  diff(Ex, 'ry');
dEx_rz =  diff(Ex, 'rz');

% dEx_vx =  diff(Ex, 'vx'); 
% dEx_vy =  diff(Ex, 'vy');
% dEx_vz =  diff(Ex, 'vz');

dEx_rx_ = subs(dEx_rx);
dEx_ry_ = subs(dEx_ry);
dEx_rz_ = subs(dEx_rz);

% dEx_vx_ = subs(dEx_vx);
% dEx_vy_ = subs(dEx_vy);
% dEx_vz_ = subs(dEx_vz);

% [23.04.2016]. Blasfemia!!
% dEx = dEx_rx_ * despl(1) + dEx_ry_ * despl(2) + dEx_rz_ * despl(3) + dEx_vx_ * ac_emis(1) + dEx_vy_ * ac_emis(2) + dEx_vz_ * ac_emis(3);
dEx = dEx_rx_ * despl(1) + dEx_ry_ * despl(2) + dEx_rz_ * despl(3);


dEy_rx =  diff(Ey, 'rx'); 
dEy_ry =  diff(Ey, 'ry');
dEy_rz =  diff(Ey, 'rz');

% dEy_vx =  diff(Ey, 'vx'); 
% dEy_vy =  diff(Ey, 'vy');
% dEy_vz =  diff(Ey, 'vz');

dEy_rx_ = subs(dEy_rx);
dEy_ry_ = subs(dEy_ry);
dEy_rz_ = subs(dEy_rz);

% dEy_vx_ = subs(dEy_vx);
% dEy_vy_ = subs(dEy_vy);
% dEy_vz_ = subs(dEy_vz);

% [23.04.2016]. Blasfemia!!
% dEy = dEy_rx_ * despl(1) + dEy_ry_ * despl(2) + dEy_rz_ * despl(3) + dEy_vx_ * ac_emis(1) + dEy_vy_ * ac_emis(2) + dEy_vz_ * ac_emis(3);
dEy = dEy_rx_ * despl(1) + dEy_ry_ * despl(2) + dEy_rz_ * despl(3);



dEz_rx =  diff(Ez, 'rx'); 
dEz_ry =  diff(Ez, 'ry');
dEz_rz =  diff(Ez, 'rz');

% dEz_vx =  diff(Ez, 'vx'); 
% dEz_vy =  diff(Ez, 'vy');
% dEz_vz =  diff(Ez, 'vz');

dEz_rx_ = subs(dEz_rx);
dEz_ry_ = subs(dEz_ry);
dEz_rz_ = subs(dEz_rz);

% dEz_vx_ = subs(dEz_vx);
% dEz_vy_ = subs(dEz_vy);
% dEz_vz_ = subs(dEz_vz);

% [23.04.2016]. Blasfemia!!
% dEz = dEz_rx_ * despl(1) + dEz_ry_ * despl(2) + dEz_rz_ * despl(3) + dEz_vx_ * ac_emis(1) + dEz_vy_ * ac_emis(2) + dEz_vz_ * ac_emis(3);
dEz = dEz_rx_ * despl(1) + dEz_ry_ * despl(2) + dEz_rz_ * despl(3);


dBx_rx =  diff(Bx, 'rx'); 
dBx_ry =  diff(Bx, 'ry');
dBx_rz =  diff(Bx, 'rz');

% dBx_vx =  diff(Bx, 'vx'); 
% dBx_vy =  diff(Bx, 'vy');
% dBx_vz =  diff(Bx, 'vz');

dBx_rx_ = subs(dBx_rx);
dBx_ry_ = subs(dBx_ry);
dBx_rz_ = subs(dBx_rz);

% dBx_vx_ = subs(dBx_vx);
% dBx_vy_ = subs(dBx_vy);
% dBx_vz_ = subs(dBx_vz);

% [23.04.2016]. Blasfemia!!
% dBx = dBx_rx_ * despl(1) + dBx_ry_ * despl(2) + dBx_rz_ * despl(3) + dBx_vx_ * ac_emis(1) + dBx_vy_ * ac_emis(2) + dBx_vz_ * ac_emis(3);
dBx = dBx_rx_ * despl(1) + dBx_ry_ * despl(2) + dBx_rz_ * despl(3);


dBy_rx =  diff(By, 'rx'); 
dBy_ry =  diff(By, 'ry');
dBy_rz =  diff(By, 'rz');

% dBy_vx =  diff(By, 'vx'); 
% dBy_vy =  diff(By, 'vy');
% dBy_vz =  diff(By, 'vz');

dBy_rx_ = subs(dBy_rx);
dBy_ry_ = subs(dBy_ry);
dBy_rz_ = subs(dBy_rz);

% dBy_vx_ = subs(dBy_vx);
% dBy_vy_ = subs(dBy_vy);
% dBy_vz_ = subs(dBy_vz);

% [23.04.2016]. Blasfemia!!
% dBy = dBy_rx_ * despl(1) + dBy_ry_ * despl(2) + dBy_rz_ * despl(3) + dBy_vx_ * ac_emis(1) + dBy_vy_ * ac_emis(2) + dBy_vz_ * ac_emis(3);
dBy = dBy_rx_ * despl(1) + dBy_ry_ * despl(2) + dBy_rz_ * despl(3);



dBz_rx =  diff(Ez, 'rx'); 
dBz_ry =  diff(Ez, 'ry');
dBz_rz =  diff(Ez, 'rz');

% dBz_vx =  diff(Ez, 'vx'); 
% dBz_vy =  diff(Ez, 'vy');
% dBz_vz =  diff(Ez, 'vz');

dBz_rx_ = subs(dEz_rx);
dBz_ry_ = subs(dEz_ry);
dBz_rz_ = subs(dEz_rz);

% dBz_vx_ = subs(dEz_vx);
% dBz_vy_ = subs(dEz_vy);
% dBz_vz_ = subs(dEz_vz);

% [23.04.2016]. Blasfemia!!
% dBz = dBz_rx_ * despl(1) + dBz_ry_ * despl(2) + dBz_rz_ * despl(3) + dBz_vx_ * ac_emis(1) + dBz_vy_ * ac_emis(2) + dBz_vz_ * ac_emis(3);
dBz = dBz_rx_ * despl(1) + dBz_ry_ * despl(2) + dBz_rz_ * despl(3);



dE = [dEx, dEy, dEz];

dB = [dBx, dBy, dBz];


% Faltaria aplicar el doppler y sacar la derivada del doppler.
doppler_invers = norm(R_vect) / (norm(R_vect) - dot(R_vect, v_vect));

dop_invers_test = subs(Dop_inv);

[E_rad, B_rad] = CampoRadiadoCargaAcelerada(R_vect, v_ini, ac_emis);


dF = dE + cross( v_vect, dB) + cross( dv_recept , B_rad);


F = E_rad + cross( v_vect, B_rad);

% Compruebo que la formula empleada sera correcta
Ex_ = subs(Ex);
Ey_ = subs(Ey);
Ez_ = subs(Ez);
E_rad_formula = [Ex_, Ey_, Ez_];

Bx_ = subs(Bx);
By_ = subs(By);
Bz_ = subs(Bz);
B_rad_formula = [Bx_, By_, Bz_];


% dDoppler_inv incluye variacion por cambio de velocidad y de posicion del
% receptor
dp = dF * doppler_invers + F * dDoppler_inv_;

% Correccion del 30.04 para incluir el efecto de la variacion del doppler
% inverso
%dE = dot(dE, v_vect) + dot(E_rad, dv_recept);
dE = dot(dE, v_vect) + dot(E_rad, dv_recept) + dot(E_rad, v_vect) * dDoppler_inv_;


Dif_Pot_EM = [dE, dp(1), dp(2), dp(3)];