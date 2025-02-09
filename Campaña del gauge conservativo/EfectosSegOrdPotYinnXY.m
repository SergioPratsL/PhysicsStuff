function ValSciGod = EfectosSegOrdPotYinnXY( Rx, Ry, v, dv_e, dv_r, despl )
% Package: Prueba de la aceleracion dual. Ultimo enfoque (esta vez si)

% Esta funcion devuelve la variacion de potencial de segudo orden debida a
% la aceleracion del emisor y del receptor.
% Estos efectos de segundo orden me han llevado por la calle de la amargura
% durante casi 3 meses pero ahora los he entendido, aunque el riesgo de
% dejarme algun termino hace que esto pueda ser una ruina

% Lo que quiero obtener (es bien dificil), es la difere entre el potencial
% que se ve tras aplicar el efecto del emisor y del receptor menos el
% potencial que se ve aplicando el efecto del emisor menos el potencial que
% se  ve aplicando solo el efecto del receptor mas el potencial que se ve
% sin aplicar ningun efecto


% dv_e es la variacion de velocidad del emisor debido a su "acelerin"
% dv_r es la variacion de velocidad del receptor 
% Tanto dv_e como dv_r estan medidas desde el SRI del emisor de forma que
% a dv_r se le ha aplicado la transformacion de velocidades mientras que a

% Ambas son infinitesimales

% v es la velocidad a la que se mueve el receptor, dado que dv_e es
% velocidad propia, se le ha de aplicar transformacion de velocidades para
% sacar la variacion sobre la velocidad v.

% dv_r es la variacion de poscion debida a la variacion de velocidad del
% receptor y al instante en que acelero

% Todo normalizado respecto a la velocidad de la luz


R = sqrt(Rx^2 + Ry^2);
Sigma = 1 / sqrt(1-v^2);

[d1_V, d1_Ax, d1_Ay, d1_Az] = PrimDerivParc(Rx, Ry, v);

Vect = [Rx, Ry, 0];
v_vect = [v, 0, 0];
[V_yinn, A_yinn] = PotencialDeYinn(Vect, v_vect);


% Primero hay que sacar los desplazamientos, rotaciones y cambios de la velocidad de cada uno 
% SIN efectos de orden 2
% R + dR_e, R + dR_r

Rot_e = dv_e(2) / (v*Sigma);
Rot_r = - dv_r(2) / v;

aux_rot_e = Rot_e * [-Ry, Rx, 0];
dR_e = - R * dv_e + aux_rot_e;

aux_rot_r = Rot_r * [-Ry, Rx, 0];
dR_r = despl + aux_rot_r;


dR_t = dR_e + dR_r;


%Ahora saco la variacion del modulo de la velocidad por efecto de ambos
%(sin orden 2)
dvel_e = - dv_e(1) / Sigma^2;

dvel_r = dv_r(1);

dvel_t = dvel_e + dvel_r;


% Para las rotaciones y boost inverso hay que saber todos los efectos de
% primer orden, pasare de los scripts anteriores porque tienen otro enfoque
dV_e_lv1 = dR_e(1) * d1_V(1) + dR_e(2) * d1_V(2) + dv_e(1) * A_yinn(1) + dv_e(2) * A_yinn(2);

dV_r_lv1 = dR_r(1) * d1_V(1) + dR_r(2) * d1_V(2);


dAx_e_lv1 = dR_e(1) * d1_Ax(1) + dR_e(2) * d1_Ax(2) - dvel_e * A_yinn(1) / v + Rot_e * A_yinn(2) + dv_e(1) * V_yinn;

dAx_r_lv1 = dR_r(1) * d1_Ax(1) + dR_r(2) * d1_Ax(2) - dvel_r * A_yinn(1) / v + Rot_r * A_yinn(2);


dAy_e_lv1 = dR_e(1) * d1_Ay(1) + dR_e(2) * d1_Ay(2) - dvel_e * A_yinn(2) / v - Rot_e * A_yinn(1) + dv_e(2) * V_yinn;

dAy_r_lv1 = dR_r(1) * d1_Ay(1) + dR_r(2) * d1_Ay(2) - dvel_r * A_yinn(2) / v - Rot_r * A_yinn(1);





% Sacar la segunda derivada de la posicion:
d2_pos_V = 3 * dot( Vect, dR_t)^2 / R^5 - 3 * dot( Vect, dR_e)^2 / R^5 - 3 * dot( Vect, dR_r)^2 / R^5;

d2_pos_Ax = d2V_pos / v;

d2_pos_Ay = 0;   % De momento no me interesa



% Sacar la segunda derivada del modulo de la velocidad
d2_vel_V = 0;        % No depende de la velocidda

d2_vel_Ax = (2/(R*v) * ( dvel_t^2 - dvel_e^2 - dvel_r^2 );

d2_vel_Ay = 0;   % De momento no me interesa



% Sacar los efectos de segundo orden de la desrotacion
d2_des_rot_Ax = Rot_e * dAy_r_lv1 + Rot_r * dAy_e_lv1;

%d_des_rot_Ay = - Rot_e * dAx_r_lv1 - Rot_r * dAx_e_lv1;
d2_des_rot_Ay = 0;       % Pq de momento paso de la Y...



% Efecto de segundo orden del boost inverso
d2_boost_inv_V = dv_e(1) * dAx_r_lv1 + dv_e(2) * dAy_r_lv1;

d2_boost_inv_Ax = dv_e(1) * dV_r_lv1;

% d2_boost_inv_Ay = dv_e(2) * dV_r_lv1;
d2_boost_inv_Ay = 0;   % Pq de momento paso de la Y... 




% Ahora los terminos directos de segudo orden, el arma definitiva del
% Vacio...

% Variacion de velocidad por modificarse el v y por tanto Sigma:
% NOTA: no hay que aplicar Sigma^2 porque dv_r(1) ya lo lleva 
dif_2_Sigma_X = 2 * v * dv_e(1) * dv_r(1);






