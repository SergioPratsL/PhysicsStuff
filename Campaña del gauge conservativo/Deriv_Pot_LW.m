function dPot_dv_fut_LW = Deriv_Pot_LW(Rfx, Rfy, a)

% Prueba piloto basada en el potencial de LW para comprobar que mi planteamiento 
% sobre como obtener el campo radiado es correcto y asi apuntalar la prueba
% del campo radiado.

% v no es necesaria pues el potencial no depende de la velocidad del emisor
% y la prueba se hace en el SRI del emisor

dPot_dv_fut_LW = [0,0,0,0];

Rf = [Rfx, Rfy, 0];

dist_f = norm(Rf);

Rf_norm = Rf / dist_f;

dPot_dv_fut_LW(1) = (1/dist_f) * dot(a,Rf_norm);

dPot_dv_fut_LW(2) = a(1) / dist_f;

dPot_dv_fut_LW(3) = a(2) / dist_f;

dPot_dv_fut_LW(4) = a(3) / dist_f;



