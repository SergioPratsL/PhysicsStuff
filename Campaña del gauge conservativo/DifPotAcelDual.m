function [Dif_Pot_A, Dif_Pot_B, Dif_Pot_Tot, Dif_Pot_Tot_en_B] = DifPotAcelDual(Rx, Ry, v, t0_b, a_A, a_B)

% Package: Prueba de la aceleracion dual

% Script final(*1) del potencial conservativo

% Devuelve las diferencias que son segundas derivadas del potencial causado
% por una aceleracion instantanea (e infinitesimal) de ambos.

% Dif_Pot_A es la variacion del potencial causada en B al recibir la luz
% que llega en tf1, relacionada con la aceleracion de A.
% Dif_Pot_B es la variacion del potencial causada en A al recibir la luz
% que llega en tf1, relacionada con la aceleracion de B.


% Rx, y Ry dan la posicion actual

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


% Efecto del campo que "A" radia en 0 y B recibe en tf1

% Variacion de velocidad de B vista desde el SRI de A:
% OJO, diria que la aceleracion debe llevar un factor 1/Sigma, pero de
% hecho no nos importa
a_B_en_A = [a_B(1)/Sigma^2, a_B(2)/Sigma, a_B(3)/Sigma];

despl_B = a_B_en_A * (tf1 - t0_b);

%Dif_Pot_B_aux = SegundDerivPotYinn_v2(Rf1(1), Rf1(2), v, a_A, despl_B);
%Dif_Pot_B_aux2 = SegundDerivPotYinn_vel(Rf1(1), Rf1(2), v, a_A, a_B_en_A);

Dif_Pot_B_deriv = EfectosSegOrdPotYinnXY(Rf1(1), Rf1(2), v, a_A, despl_B, a_B_en_A);


% Sacar la diferencia de la fuerza hecha por campo radiado
[dE_B, dB_B] = DerivCampRad( Rf1(1), Rf1(2), a_A, despl_B );

[E_B, B_B] = CampoRadiadoCargaAcelerada(Rf1, v_0, a_A);

F_B = E_B + cross(v_vect, B_B);

% La variacion de la fuerza depende de la variacion de E y B pero tambien
% de la variacion de la velocidad del receptor
dF_B = dE_B + cross( v_vect, dB_B ) + cross(a_B_en_A, B_B);

% La variacion del trabajo hecho dependera no solo de la variacion de la
% fuerza sino tambien de la de la  velocidad del receptor
dPotencia_B = dot(dF_B, v_vect) + dot(F_B, a_B_en_A);
%dPotencia_B = dot(dE_B, v_vect) + dot(F_B, a_B_en_A);

% 04.01.2016. El puto Doppler :D
Dopler_B = 1 * (1 - v * Rf1(1) / tf1);

% 03.01.2016. Cuestiono el signo menos
%Dif_Pot_B_EM = [-dPotencia_B, -dF_B(1), -dF_B(2), -dF_B(3)];
% Y añado el Doppler
Dif_Pot_B_EM = [dPotencia_B, dF_B(1), dF_B(2), dF_B(3)] / Dopler_B

Dif_Pot_B = Dif_Pot_B_deriv + Dif_Pot_B_EM;



% Efecto del campo que "A" radia en 0 y B recibe en tf1

% Para el efecto en tf2 es mas chungo... primero hay que obtener las
% aceleraciones en el SRI de B lo cual ya tiene tela
a_A_en_B = [a_A(1)/Sigma^2, a_A(2)/Sigma, a_A(3)/Sigma];

% Tb hay que calcular el tiempo que ha estado viajando, que es entre tf2 y
% el instante 0 evaluado desde B, que ya no es cero.
% ¿Cuanto tiempo ha pasado en el SRI de B desde que A acelero hasta que A
% recibe la luz de tf2?? --> Dos sucesos sobre "A" --> mirar de dilatar el
% tiempo
% La luz de t0_b llega en t0_b mas norm(R2), por tanto lo  tenemos
intervalo = Sigma * tf2;

despl_A = a_A_en_B * intervalo;

%Dif_Pot_A_en_B_aux = SegundDerivPotYinn_v2(Rf2_b(2), Rf2_b(3), -v, a_B, despl_A);
%Dif_Pot_A_en_B_aux2 = SegundDerivPotYinn_vel(Rf2_b(2), Rf2_b(3), -v, a_B, a_A_en_B);

Dif_Pot_A_en_B_deriv = EfectosSegOrdPotYinnXY(Rf2_b(2), Rf2_b(3), -v, a_B, despl_A, a_A_en_B);

% Sacar la diferencia de la fuerza hecha por campo radiado
[dE_A, dB_A] = DerivCampRad( Rf2_b(2), Rf2_b(3), a_B, despl_A );

Rf2_aux = [Rf2_b(2), Rf2_b(3), Rf2_b(4)];

[E_A, B_A] = CampoRadiadoCargaAcelerada(Rf2_aux, v_0, a_B);

F_A = E_A + cross(-v_vect, B_A);

% La variacion de la fuerza depende de la variacion de E y B pero tambien
% de la variacion de la velocidad del receptor
dF_A = dE_A + cross( -v_vect, dB_A ) + cross(a_A_en_B, B_A);

% La variacion del trabajo hecho dependera no solo de la variacion de la
% fuerza sino tambien de la de la  velocidad del receptor
dPotencia_A = dot(dF_A, -v_vect) + dot(F_A, a_A_en_B);
%dPotencia_A = dot(dE_A, -v_vect) + dot(F_A, a_A_en_B);

% 04.01.2016. El puto Doppler :D
dist_A = norm(Rf2_aux);
Dopler_A = 1 * (1 - (-v) * Rf2_b(2) / dist_A);

% 03.01.2016. Cuestiono el signo menos
%Dif_Pot_A_EM_en_B = [-dPotencia_A, -dF_A(1), -dF_A(2), -dF_A(3)];
Dif_Pot_A_en_B_EM = [dPotencia_A, dF_A(1), dF_A(2), dF_A(3)] / Dopler_A;

Dif_Pot_A_en_B = Dif_Pot_A_en_B_deriv + Dif_Pot_A_en_B_EM

% Sumar los dos factores para sacar la diferencia
Dif_Pot_A  = Boost( Dif_Pot_A_en_B, -v_vect );


Dif_Pot_Tot = Dif_Pot_A + Dif_Pot_B;

% Ver el valor desde B
Dif_Pot_Tot_en_B = Boost( Dif_Pot_Tot, v_vect );


auxy = 1;

























% (*1)... Nunca se sabe...

