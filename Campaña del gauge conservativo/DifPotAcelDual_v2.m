function [Dif_Pot_A, Dif_Pot_B, Dif_Pot_Tot, Dif_Pot_Tot_en_B] = DifPotAcelDual_v2(Rx, Ry, v, t0_b)

% Package: Prueba de la aceleracion dual

% Script final(*1) del potencial conservativo

% Devuelve las diferencias que son segundas derivadas del potencial causado
% por una aceleracion instantanea (e infinitesimal) de ambos.

% Dif_Pot_A es la variacion del potencial causada en B al recibir la luz
% que llega en tf1, relacionada con la aceleracion de A.
% Dif_Pot_B es la variacion del potencial causada en A al recibir la luz
% que llega en tf1, relacionada con la aceleracion de B.


% Rx, y Ry dan la posicion actual ¡de B respecto a A!

% v la velocidad

% t_0b indica el instante en que acelera B, siendo que A acelera en el
% instante 0

% a_A es la aceleracion vectorial de A vista desde el SRI de A
% a_B es la aceleracoin de B vista desde el SRI de B
% espero que Thomas no lo joda todo...

% Todo normalizado respecto la velocidad de la luz



% Primero, obtener posiciones retardadas


R1 = [Rx, Ry, 0]; 
v_vect = [v,0,0];

Sigma = 1 / sqrt(1 - v^2);

Rr1 = PostRetFut( R1, v, 'r');
Rr1= - Rr1;
tr1 = - norm(Rr1);

Rf1 = PostRetFut( R1, v, 'f');
tf1 = norm(Rf1);

% Validamos que t0_b este en el intervalo enrte 
if t0_b <= tr1 || t0_b >= tf1
    disp('Tiempo fuera de intervalo, cancelar')
    return
end

v_0 = [0, 0, 0];

% Obtener la posicion futura de "B" cuando se recibe la luz que envia en
% t0_b. Los signos son para hacer el vector de "B" a "A"
Rx2 = - Rx - t0_b * v;
Ry2 = - Ry;
R2 = [Rx2, Ry2, 0];
t2 = norm(R2);

tf2 = t0_b + t2; 

% Al estar en el SRI en que no se mueve "A", la distancia retardada y la futura solo cambian de signo, 
% el tiempo en ambas siempre es positivo, pero el tiempo afectada de manera
% opuesta a un lado y al otro.


R2_t = [t2, Rx2, Ry2, 0];


% Boost al sistema de "B" para ver la posicion futura de "A" vista desde el
% SRI de B
Rf2_b = Boost( R2_t,  v_vect );


% Efecto del campo inducido que provoca la fuerz

% Hemos de calcular a_A como la fuerza vista desde el SRI del receptor, es
% mas facil aunque no puedas usar Coulomb, se ha de tener en cuenta la
% posicion pasada
[E_a, B_a] = CampoInducido_sin_unidades(Rr1, v_vect);

% Fuerza vista desde el SRI propio. Desde el SRI receptor no usar -v_vect sino v_0!!
a_A = FLorentz(v_0, E_a, B_a);

a_A_en_B = [a_A(1)/Sigma^2, a_A(2)/Sigma, a_A(3)/Sigma];

% Para obtener A_B  necesitamos la posicion pasada que ve "B"...
% Primero lo sacamos en el SRI de A (vA = 0)
Rrx2 = (Rx + t0_b * v);
Rry2 = Ry;
% Rr2_en_A = [Rx2, Ry2, 0];
Rr2_en_A = [Rrx2, Rry2, 0];
tr2_en_A = norm(Rr2_en_A);

Rr2_t_en_A = [tr2_en_A, Rr2_en_A(1), Rr2_en_A(2), Rr2_en_A(3)];

Rr2_t = Boost( Rr2_t_en_A, v_vect );
Rr2 = [Rr2_t(2), Rr2_t(3), Rr2_t(4)];

[E_b, B_b] = CampoInducido_sin_unidades(Rr2, -v_vect);

% Desde el SRI receptor no usar -v_vect sino v_0!!
a_B = FLorentz(v_0, E_b, B_b);


% Variacion de velocidad de B vista desde el SRI de A:
% OJO, diria que la aceleracion debe llevar un factor 1/Sigma, pero de
% hecho no nos importa pues solo importa para la radiacion pero esa siempre
% la sacamos en el SRI propio.
a_B_en_A = [a_B(1)/Sigma^2, a_B(2)/Sigma, a_B(3)/Sigma];

R_aux = [Rf1(1), Rf1(2), Rf1(3)];
R_aux = R_aux / norm(R_aux);
despl_B = a_B_en_A * (tf1 - t0_b) + (tf1 - t0_b) * dot(a_B_en_A, R_aux)/(1-dot(R_aux, v_vect)) * v_vect;

Dif_intervalo_A = (tf1 - t0_b) * dot(a_B_en_A, R_aux)/(1-dot(R_aux, v_vect));

%Dif_Pot_B_aux = SegundDerivPotYinn_v2(Rf1(1), Rf1(2), v, a_A, despl_B);
%Dif_Pot_B_aux2 = SegundDerivPotYinn_vel(Rf1(1), Rf1(2), v, a_A, a_B_en_A);

% A_B ya no se mide desde el emisor (ya no es a_B_en_A)
Dif_Pot_B_deriv = efectos_segundo_orden_Scigod(Rf1, v_vect, a_A, a_B, despl_B, Dif_intervalo_A);

Dif_Pot_B_EM = EfectosCamposRadiado(Rf1, v_vect, despl_B, a_A, a_B_en_A);

Dif_Pot_B = Dif_Pot_B_deriv + Dif_Pot_B_EM 


% Efecto del campo que "A" radia en 0 y B recibe en tf1

% Para el efecto en tf2 es mas chungo... primero hay que obtener las
% aceleraciones en el SRI de B lo cual ya tiene tela

% Tb hay que calcular el tiempo que ha estado viajando, que es entre tf2 y
% el instante 0 evaluado desde B, que ya no es cero.
% ¿Cuanto tiempo ha pasado en el SRI de B desde que A acelero hasta que A
% recibe la luz de tf2?? --> Dos sucesos sobre "A" --> mirar de dilatar el
% tiempo
% La luz de t0_b llega en t0_b mas norm(R2), por tanto lo  tenemos
intervalo = Sigma * tf2;

Rf2 = [Rf2_b(2), Rf2_b(3), Rf2_b(4)];
R_aux2 = Rf2 / norm(Rf2);

despl_A = a_A_en_B * intervalo + intervalo * dot(a_A_en_B, R_aux2)/(1-dot(R_aux2, -v_vect)) * (-v_vect);

Dif_intervalo_B = (tf1 - t0_b) * dot(a_A_en_B, R_aux2)/(1-dot(R_aux2, -v_vect));

%Dif_Pot_A_en_B_aux = SegundDerivPotYinn_v2(Rf2_b(2), Rf2_b(3), -v, a_B, despl_A);
%Dif_Pot_A_en_B_aux2 = SegundDerivPotYinn_vel(Rf2_b(2), Rf2_b(3), -v, a_B, a_A_en_B);

% Dif_Pot_A_en_B_deriv = SegundDerivPotYinnUnif(Rf2_b(2), Rf2_b(3), -v, a_B, a_A_en_B, despl_A);
Dif_Pot_A_en_B_deriv = efectos_segundo_orden_Scigod(Rf2, -v_vect, a_B, a_A, despl_A, Dif_intervalo_B);

Dif_Pot_A_en_B_EM = EfectosCamposRadiado(Rf2, -v_vect, despl_A, a_B, a_A_en_B);

Dif_Pot_A_en_B = Dif_Pot_A_en_B_deriv + Dif_Pot_A_en_B_EM 

% Sumar los dos factores para sacar la diferencia
Dif_Pot_A  = Boost( Dif_Pot_A_en_B, -v_vect );




Dif_Pot_Tot = Dif_Pot_A + Dif_Pot_B;

% Ver el valor desde B
Dif_Pot_Tot_en_B = Boost( Dif_Pot_Tot, v_vect );


auxy = 1;


