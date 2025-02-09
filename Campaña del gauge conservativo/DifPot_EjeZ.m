
function [Pot_A, Pot_B, F_rad, Dif] = DifPot_EjeZ( Rx, Ry,  v)

% Rx y Ry representan la posicion actual, esto es asi porque hay que
% trabajar tanto con la distancia retardada como con la futura.
% Rx, Ry van desde el emisor (que acelera) hasta el receptor

% La distancia retardada para el efecto del potencial del que acelera "A"
% la distancia futura para el efecto del potencial del que no acelera "B"
% asi como para el campo radiado

% v es un numero que representa la velocidad que de "B" sobre el eje X
% normalizada a la velocidad de la luz. Estamos en SRI de "A"

% Pot_A: variacion del potencial sobre el que acelera
% Pot_B: variacion del potencial sobre el que no acelera
% F_rad: potencia y fuerza del campo radiado sobre el que no acelera
% Dif: aunque pueda engañar el nombre, es la suma de los tres anteriores

% El acelerin en direccion Y, dv, se normaliza a 1... y se asume en sentido
% positivo...


R = [Rx, Ry, 0]; 
v_vect = [v,0,0];

Sigma = 1 / sqrt(1 - v^2);

Rr = PostRetFut( R, v, 'r');
% Hay que acerlo asi para que vaya bienm si inviertes en el paso anterior
% la lias y obtienes el vector futuro pero invertido.
Rr = - Rr;

Rf = PostRetFut( R, v, 'f');


% Hago Boost del vector Rr para evaluar el potencial en el SRI de "B"
tr = norm(Rr);
Rr_t = [tr, Rr];        % '_t' porque incluye el tiempo

Rr_B_t = Boost(Rr_t, v_vect);


% La rotacion tendra el factor 1/Sigma al estar en B. Se rota entre el eje
% X y el Z puesto que la velocidad era en X...
dRot_1 = 1 / (v*Sigma);

dPot_A_dx = DerivadasPotencialYinn( Rr_B_t(2), Rr_B_t(3), -v, 'x');

dPot_A_dz = DerivadasPotencialYinn( Rr_B_t(2), Rr_B_t(3), -v, 'z');



% Variacion debido a como la primera rotacion altera Rz (Rx no al ser Rz =
% 0)
dV_dRot_A = dRot_1 * Rr_B_t(2) * dPot_A_dz(1);

dAx_dRot_A = dRot_1 * Rr_B_t(2) * dPot_A_dz(2);

dAy_dRot_A = dRot_1 * Rr_B_t(2) * dPot_A_dz(3);

dAz_dRot_A = dRot_1 * Rr_B_t(2) * dPot_A_dz(4);

d4V_dRot_A = [dV_dRot_A, dAx_dRot_A, dAy_dRot_A, dAz_dRot_A];


Rr_en_B = [Rr_B_t(2), Rr_B_t(3), Rr_B_t(4)];

[V_A_desde_B, A_A_desde_B] = PotencialDeYinn(Rr_en_B, -v_vect);

% Variacion debido a la segunda rotacion para volver al plano original
dV_desrot_A = 0;

dAx_desrot_A = dRot_1 * A_A_desde_B(3);

dAy_desrot_A = 0;

dAz_desrot_A = - dRot_1 * A_A_desde_B(1);

d4V_desrot_A = [dV_desrot_A, dAx_desrot_A, dAy_desrot_A, dAz_desrot_A];


dPot_A_desde_B = d4V_dRot_A + d4V_desrot_A;

% A esta  variacion de potencial hay que hacerle un boost para que vuelva
% al SRI origen

Pot_A = boost(dPot_A_desde_B, -v_vect);




% Obtener la variacion de potencial de A sobre B

% Primero  se hace un miniboost que modifica Rz en  - dv* R, la
% modificacion del tiempo no afecta porque el potencial depende solo de la
% distancia
dist_f = norm(Rf);

dPot_B_dx = DerivadasPotencialYinn( Rf(1), Rf(2), v, 'x');

dPot_B_dz = DerivadasPotencialYinn( Rf(1), Rf(2), v, 'z');

dPot_B_mini_boost = - dist_f * dPot_B_dz;

% Al no ser en direccion X la velocidad vista desde A', hay que hacer la
% mini rotacion En este caso la rotacion no llevara el factor 1/Sigma 
dRot_2 = 1 / v;

dV_dRot_B = dRot_2 * Rf(1) * dPot_B_dz(1);

dAx_dRot_B = dRot_2 * Rf(1) * dPot_B_dz(2);

dAy_dRot_B = dRot_2 * Rf(1) * dPot_B_dz(3);

dAz_dRot_B = dRot_2 * Rf(1) * dPot_B_dz(4);

d4V_dRot_B = [dV_dRot_B, dAx_dRot_B, dAy_dRot_B, dAz_dRot_B];


% Desrotamos sobre el potencial
[V_A_desde_A, A_A_desde_A] = PotencialDeYinn(Rf, v_vect);

dV_desrot_B = 0;

dAx_desrot_B = dRot_2 * A_A_desde_A(3);

dAy_desrot_B = 0;

dAz_desrot_B = - dRot_2 * A_A_desde_A(1);

d4V_desrot_B = [dV_desrot_B, dAx_desrot_B, dAy_desrot_B, dAz_desrot_B];


% Deshacemos el mini-boost
dV_boost_inv_B = A_A_desde_A(3);  % Entiendo que sera 0...

dAx_boost_inv_B = 0;

dAy_boost_inv_B = 0;

dAz_boost_inv_B = V_A_desde_A;

dPot_B_boost_inv = [dV_boost_inv_B, dAx_boost_inv_B, dAy_boost_inv_B, dAz_boost_inv_B];


 % El potencial es la suma de todo lo anterior:
Pot_B = dPot_B_mini_boost + d4V_dRot_B + d4V_desrot_B + dPot_B_boost_inv;



% Campo radiado:
v_0 = [0, 0, 0];    % Se evalua desde el emisor.
dv_aux = [0, 0, 1]; 

[E, B] = CampoRadiadoCargaAcelerada(Rf, v_0, dv_aux);

Dopler = 1 / (1 - v * Rf(1) / dist_f);

F = E + cross(v_vect, B);

F_rad = [dot(F, v_vect), F] * Dopler;



Dif = Pot_A + Pot_B + F_rad;

