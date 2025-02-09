
function [Pot_A, Pot_B, F_rad, Dif] = DifPot_EjeX( Rx, Ry,  v)

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

% El acelerin en direccion X, dv, se normaliza a 1... y se asume en sentido
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

% Obtener las derivadas segun las formulas que tengo apuntadas
dV_dv = 0;

% Ojo, "B" ve a "A" con velocidad "-v"

dAx_dv = - 1 / ((-v)^2 * Rr_B_t(1));    % Rr_b_t(1) es el tiempo, es decir, la distancia

dAy_dv = Rr_B_t(2) / ((-v)^2 * Rr_B_t(1) * Rr_B_t(3));

% Hay que incluir la derivada de la constante en dDy (es constante respecto
% al espacio pero no respecto a la velocidad).
dC_dv =  abs(-v) / ((-v)^3 * Rr_B_t(3));
dAy_dv = dAy_dv + dC_dv;

dAz_dv = 0;

% La aceleracion en este entorno es menor en un factur 1/Sigma^2
dPot_A_en_B = [dV_dv, dAx_dv, dAy_dv, dAz_dv] / Sigma^2;

% PRIMERA DE LAS VARIABLES, LO QUE CAMBIA EL POTENCIAL DE A P(CREADO POR B)
% POR LA ACELERACION DE A
Pot_A = Boost( dPot_A_en_B,  -v_vect );


% Un boost infinitesimal implica desplazar -dv*R en direccion X con lo que
% neesitamos la derivada en direccion X, tambien hay derivada de la
% velocidad y luego un boost de -dv sobre el vector 4-potencial
% La derivada temporal es cero asi que no nos preocupa.

dPot_B_dx = DerivadasPotencialYinn( Rf(1), Rf(2), v, 'x');

dist_f = norm(Rf);

dV_dv = 0;
dAx_dv = - 1 / ( v^2 * dist_f );
dAy_dv = Rf(1) / ( v^2 * dist_f * Rf(2) );

dC_dv = abs(v) / ( v^3 * Rf(2));
dAy_dv = dAy_dv + dC_dv;

dAz_dv = 0;

% Si acelero el origen la aceleracion es en sentido contrario (signo
% menos)! Tambien debe aplicars el factor 1/Sigma^2
dPot_B_dv = - [dV_dv, dAx_dv, dAy_dv, dAz_dv] / Sigma^2;


[V_B, A_B] = PotencialDeYinn(Rf, v_vect);

dPot_Boost_final =  [A_B(1), V_B, 0, 0];

% SEGUNDA DE LAS VARIABLES: LO QUE CAMBIA EL POTENCIAL DE B (CREADO POR A)
% POR LA ACELERACION DE A
Pot_B = - dist_f * dPot_B_dx + dPot_B_dv + dPot_Boost_final;

% Prueba porque no se si los boosts aplican.
% Pot_B = dPot_B_dv;

% Campo radiado:
v_0 = [0, 0, 0];    % Se evalua desde el emisor.
dv_aux = [1, 0, 0]; 

[E, B] = CampoRadiadoCargaAcelerada(Rf, v_0, dv_aux);

% 04.01.2016. Cambio inocuo pues la definicion del doppler es sobre
% frecuencia no sobre tiempo y por ello se han de invertir los signos!
Dopler = 1 * (1 - v * Rf(1) / dist_f);

F = E + cross(v_vect, B);

% 04.01.2016. Cambio inocuo "/ Dopler" y no "* Dopler"
F_rad = [dot(F, v_vect), F] / Dopler;




Dif = Pot_A + Pot_B + F_rad;


ZZ = 1;