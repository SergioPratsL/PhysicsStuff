function [E, B] = CampoRadiadoDerivadas_LW( Rfx, Rfy, a);


Rf = [Rfx, Rfy, 0];

dist_f = norm(Rf);

Rf_norm = Rf / dist_f;

% Estas derivadas son por unidad de  tiempo.
dPot_dv_fut_LW = Deriv_Pot_LW(Rfx, Rfy, a);

dPot_dv_fut_LW_dt = dPot_dv_fut_LW;

dPot_dv_fut_LW_dx = - Rf_norm(1)  * dPot_dv_fut_LW;

dPot_dv_fut_LW_dy = - Rf_norm(2) * dPot_dv_fut_LW;

% Dado que Rfz = 0
dPot_dv_fut_LW_dz = [0,0,0,0];


% Obtenemos el campo
Ex = - dPot_dv_fut_LW_dt(2) - dPot_dv_fut_LW_dx(1);
Ey = - dPot_dv_fut_LW_dt(3) - dPot_dv_fut_LW_dy(1);
Ez = - dPot_dv_fut_LW_dt(4) - dPot_dv_fut_LW_dz(1);

E = [Ex, Ey, Ez];


Bx = dPot_dv_fut_LW_dy(4) - dPot_dv_fut_LW_dz(3);
By = dPot_dv_fut_LW_dz(2) - dPot_dv_fut_LW_dx(4);
Bz = dPot_dv_fut_LW_dx(3) - dPot_dv_fut_LW_dy(2);

B = [Bx, By, Bz];


