function [E, B] = CampoRadiadoYinn( Rx, Ry, v, a);

% Esta funcion obtiene la distancia pasada y la futura a partir de la
% actual y luego saca el campo radiado



R = [Rx, Ry, 0]; 
v_vect = [v,0,0];

% Le pongo otro nombre en vez de Rr, no sea que el Matlab no respete las
% variables en diferentes funciones :S.
R_ret = PostRetFut( R, v, 'r');
R_ret = - R_ret;

R_fut = PostRetFut( R, v, 'f');

dist_fut = norm(R_fut);

R_fut_norm = R_fut / dist_fut;


dist_ret = norm(R_ret);

R_ret_norm = R_ret / dist_ret;


% Derivadas del potencial
dPot_dv_fut = Deriv_Pot_fut(R_fut(1), R_fut(2), v, a);


% "Campo de aceleracion"
dPot_dv_ret = Deriv_Pot_ret(R_ret(1), R_ret(2), v, a);


% Mal planteamiento!
% dPot_total = dPot_dv_fut + dPot_dv_ret;
%dPot_total = dPot_dv_fut;  % Peor aun :S

Coef_dx_f = - R_fut_norm(1);

%16.12.2015. Tienen que coincidir porque el instante en que se recibe es el
%mismo que el que se emite!!
%Coef_dx_r =  (R_ret_norm(1) - R_fut_norm(1) ) / (1 - v * R_ret_norm(1)) + R_fut_norm(1);
Coef_dx_r = Coef_dx_f;


Coef_dy_f = - R_fut_norm(2);

%16.12.2015. Tienen que coincidir porque el instante en que se recibe es el
%mismo que el que se emite!!
%Coef_dy_r = (R_ret_norm(2) - R_fut_norm(2)) / (1 - v * R_ret_norm(1)) + R_fut_norm(2);
Coef_dy_r = Coef_dy_f;


dPot_total_dt = dPot_dv_fut + dPot_dv_ret;

dPot_total_dx = Coef_dx_f * dPot_dv_fut + Coef_dx_r * dPot_dv_ret;

dPot_total_dy = Coef_dy_f * dPot_dv_fut + Coef_dy_r * dPot_dv_ret;;

% Dado que Rfz = 0
dPot_total_dz = [0,0,0,0];




% Obtenemos el campo
Ex = - dPot_total_dt(2) - dPot_total_dx(1);
Ey = - dPot_total_dt(3) - dPot_total_dy(1);
Ez = - dPot_total_dt(4) - dPot_total_dz(1);

E = [Ex, Ey, Ez];


Bx = dPot_total_dy(4) - dPot_total_dz(3);
By = dPot_total_dz(2) - dPot_total_dx(4);
Bz = dPot_total_dx(3) - dPot_total_dy(2);

B = [Bx, By, Bz];





