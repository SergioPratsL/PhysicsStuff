
function [Eo, Bo] = CampoRadiadoCargaAcelerada( R, v, dv)
% Campo, quitando un monton de unidades generado a partir de una carga
% que acelera. Se excluye el campo inducido. Dv es la variacion de la
% velocidad por unidad de tiempo, se separa de la aceleracion en cuanto
% depende del SRI
% R: posicion retardad de la carga que genera el campo en la posicion (x,y,z)
% v: veclodad de la carga (vx, vy,vz)
% dv: variacion de la velocidad por unidad de tiempo (aceleracion), medida
% desde el sistema del laboratorio (no es aceleracion propia) y descartando
% el efecto doppler espacial. Es dv/dt' (Ojo, hay un Sigma que baila!)
% velocidad normalizada (c=1)
% Carga normalizada (q=1) y constante de permitividad por 4pi normalizado
% tambien (solo se tienen en cuenta las consideraciones geometricas y de
% relatividad).

% Voy a sufrir para hacer cuadrar esta formula con el doble boost


dist = norm(R); % tambien es el tiempo (con c=1)
R_dot_v = dot(R, v);

R_actual = R - dist * v; 
dist_actual = norm(R_actual);

S = dist - R_dot_v;


Sigma = fGamma(v);
% Ojo con este Sigma, no esta en la wikipedia pero creo que debe estar
% Lo quito, ya vigilare luego como hacer la conversion
% dv_lab = dv * (1 - R_dot_v / dist );
% dv_lab = dv * Sigma * (1 - R_dot_v / dist );

% Pasando totalmente de la wikipedia y basandome en mis datos
% dv_lab = dv / (1 - R_dot_v / dist );
%dv_lab = dv / Sigma;
% Tiene mas sentido que se haga la conversion en el script principal
dv_lab = dv;

v1 = cross( R_actual, dv_lab);

v2 = cross(R, v1);

% Eo = v2 / dist_actual^3;

% 19.04.2014. Pero vaya ruina que no hubiera visto esto!
Eo = v2 / S^3;

Bo = cross(R, Eo) / dist;






