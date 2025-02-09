function Deriv = efectos_segundo_orden_Scigod( R_vect, v_vect, dv_e, dv_r, despl, dif_interval )
% Package: Prueba de la aceleracion dual

% Esta funcion devuelve la variacion de potencial de segudo orden debida
% los efectos cruzados de la aceleracion del emisor y del receptor, 
% en los momentos evaluados
% Devuelve el efecto de segundo orden provocado por el cambio de veloidad 
% y de posicion del receptor sobre el cambio de potencia que provoca detectar 
% que la velocidad del emisor ha cambiado

% dv_e es la variacion de velocidad propia del emisor debido a su "acelerin"
% dv_r es la variacion de velocidad propia del receptor 
% Tanto dv_e y dv_r se corresponden a aceleraciones propias (esto cambia
% respecto a los modelos anteriores que castigaban a dv_r).

% Ambas son infinitesimales

% v_vect es la velocidad vectorial a la que se mueve el receptor, dado que dv_e es
% velocidad propia, se le ha de aplicar transformacion de velocidades para
% sacar la variacion sobre la velocidad v.

% dv_r es la variacion de poscion debida a la variacion de velocidad del
% receptor y al instante en que acelero

% Todo normalizado respecto a la velocidad de la luz

R = norm(R_vect);
R_norm = R_vect / R;
v = norm(v_vect);

Sigma = 1 / sqrt(1-v^2);


dv_r_red = DifVelRed( dv_r, v_vect );

dv_e_red = DifVelRed( dv_e, v_vect );

% Quiza no es necesario con dv_r_red y dv_e_red...
%dv_dvr = MatrizConvierteVelocidades(v_vect);
%dv_dve = - dv_dvr;

% Segunda derivada de la velocidad segun la formula de adicion de
% velocidades
v_ini = [0, 0, 0];
[dvx_dv2, dvy_dv2, dvz_dv2] = SecDerivLeyAddVeloc(v_vect, v_ini);

% Para sacar la variacion de segundo orden, se debe hacer con las
% variaciones de velocidad del SRI en que se evalua, es decir, con dv_e
% (sin reducir) y con dv_r_red...
dif_v_x = dv_r_red * ( dvx_dv2 * dv_e');
dif_v_y = dv_r_red * ( dvy_dv2 * dv_e');
dif_v_z = dv_r_red * ( dvz_dv2 * dv_e');

dif_v_orden2 = [dif_v_x, dif_v_y, dif_v_z];

% [dv_dvrdve, dv_dvedvr] = DifVelocOrd2( v_vect, dv_e, dv_r );

[dRx_dvdr, dRy_dvdr, dRz_dvdr] = SecDeriv_R(R_vect, v_ini);

dif_R_x = despl * ( dRx_dvdr * dv_e');
dif_R_y = despl * ( dRy_dvdr * dv_e');
dif_R_z = despl * ( dRz_dvdr * dv_e');

dif_R_orden2 = [dif_R_x, dif_R_y, dif_R_z];

dif_R_x_trasp = despl * ( dRx_dvdr' * dv_e');
dif_R_y_trasp = despl * ( dRy_dvdr' * dv_e');
dif_R_z_trasp = despl * ( dRz_dvdr' * dv_e');

dif_R_orden2_trasp = [dif_R_x_trasp, dif_R_y_trasp, dif_R_z_trasp];


% Potencial de Yinn:
[V, A] = PotencialDeYinnAmpliado(R_vect, v_vect);

% Primera derivada espacial
[dV_dr, dAx_dr, dAy_dr, dAz_dr] = PrimDerivParcEnh(R_vect, v_vect); 

% Primera derivada de velocidad
[dV_dv, dAx_dv, dAy_dv, dAz_dv] = PrimeraDerivVelEnh(R_vect, v_vect);

% Segunda derivada espacial
[dV_dr2, dAx_dr2, dAy_dr2, dAz_dr2] = SecDerivParcEnh(R_vect, v_vect);

% Segunda derivada de la velocidad
[dV_dv2, dAx_dv2, dAy_dv2, dAz_dv2] = SecDerivParcVelEnh(R_vect, v_vect);

% Segunda derivada espacial y de velocidad
[dV_drdv, dAx_drdv, dAy_drdv, dAz_drdv] = SecDerivParc_Rv_Enh(R_vect, v_vect);




% Primer termino

Dif1_V = - dot(dif_R_orden2 , dV_dr)  - R * despl * (dV_dr2 * dv_e');
Dif1_Ax = - dot(dif_R_orden2 , dAx_dr) - R * despl * (dAx_dr2 * dv_e');
Dif1_Ay = - dot(dif_R_orden2 , dAy_dr) - R * despl * (dAy_dr2 * dv_e');
Dif1_Az = - dot(dif_R_orden2 , dAz_dr) - R * despl * (dAz_dr2 * dv_e');


% Segundo termino. Signo menos porque dv/dv_e tiene signo menos. Ojo!!
Dif2_V = - despl * (dV_drdv * dv_e_red');
Dif2_Ax = - despl * (dAx_drdv * dv_e_red');
Dif2_Ay = - despl * (dAy_drdv * dv_e_red');
Dif2_Az = - despl * (dAz_drdv * dv_e_red');


% Tercer termino, debe trasponerse la matriz pues dV_drdv toca primero la
% componente espacial y no la temporal (componente espacial que viene del
% boost R*dv_e'
Dif3_V =  - R * dv_r_red * (dV_drdv' * dv_e');
Dif3_Ax = - R * dv_r_red * (dAx_drdv' * dv_e');
Dif3_Ay = - R * dv_r_red * (dAy_drdv' * dv_e');
Dif3_Az = - R * dv_r_red * (dAz_drdv' * dv_e');


% Cuarto termino. Signo menos al prmimer dv/dv_e
Dif4_V = - dv_r_red * (dV_dv2 * dv_e_red') + dot(dif_v_orden2, dV_dv);
Dif4_Ax = - dv_r_red * (dAx_dv2 * dv_e_red') + dot(dif_v_orden2, dAx_dv);
Dif4_Ay = - dv_r_red * (dAy_dv2 * dv_e_red') + dot(dif_v_orden2, dAy_dv);
Dif4_Az = - dv_r_red * (dAz_dv2 * dv_e_red') + dot(dif_v_orden2, dAz_dv);


% Quinto termino
vect_aux = dv_e(1) * dAx_dr + dv_e(2) * dAy_dr + dv_e(3) * dAz_dr;
Dif5_V = dot( despl, vect_aux);


% Surge un dilema mas sobre como colocar los terminos de las matrices
vect_aux2 = despl * (dV_dr' * dv_e);
Dif5_Ax = vect_aux2(1);
Dif5_Ay = vect_aux2(2);
Dif5_Az = vect_aux2(3);


% Sexto termino (muy relacionado con el quinto)
vect_aux3 = dv_e(1) * dAx_dv + dv_e(2) * dAy_dv + dv_e(3) * dAz_dv;
Dif6_V = dot( dv_r_red, vect_aux3);

vect_aux4 = dv_r_red * (dV_dv' * dv_e);    % Dara cero
Dif6_Ax = vect_aux4(1);
Dif6_Ay = vect_aux4(2);
Dif6_Az = vect_aux4(3);


% 23.04.2016. No aplican los  terminos del 7 al 12!!

% % Septimo termino
% % aux2 = dot(dif_R_orden2_trasp , dAx_dr)
% Dif7_V =  dot(dif_R_orden2_trasp , dV_dr) - R * dv_e * (dV_dr2 * despl');
% Dif7_Ax = dot(dif_R_orden2_trasp , dAx_dr) - R * dv_e * (dAx_dr2 * despl');
% Dif7_Ay = dot(dif_R_orden2_trasp , dAy_dr) - R * dv_e * (dAy_dr2 * despl');
% Dif7_Az = dot(dif_R_orden2_trasp , dAz_dr) - R * dv_e * (dAz_dr2 * despl');
% 
% 
% % Octavo termino
% Dif8_V = - dv_e_red * (dV_drdv * despl');
% Dif8_Ax = - dv_e_red * (dAx_drdv * despl');
% Dif8_Ay = - dv_e_red * (dAy_drdv * despl');
% Dif8_Az = - dv_e_red * (dAz_drdv * despl');
% 
% 
% % Noveno termino
% Dif9_V = - R * dv_e * (dV_drdv * dv_r_red');
% Dif9_Ax = - R * dv_e * (dAx_drdv * dv_r_red');
% Dif9_Ay = - R * dv_e * (dAy_drdv * dv_r_red');
% Dif9_Az = - R * dv_e * (dAz_drdv * dv_r_red');
% 
% 
% % Decimo termino
% Dif10_V = - dv_e_red * (dV_dv2 * dv_r_red') + dot(dif_v_orden2, dV_dv);
% Dif10_Ax = - dv_e_red * (dAx_dv2 * dv_r_red') + dot(dif_v_orden2, dAx_dv);
% Dif10_Ay = - dv_e_red * (dAy_dv2 * dv_r_red') + dot(dif_v_orden2, dAy_dv);
% Dif10_Az = - dv_e_red * (dAz_dv2 * dv_r_red') + dot(dif_v_orden2, dAz_dv);
% 
% 
% % Undecimo termino
% vect_aux = despl(1) * dAx_dr + despl(2) * dAy_dr + despl(3) * dAz_dr;
% Dif11_V = dot( dv_e, vect_aux);
% 
% % Surge un dilema mas sobre como colocar los terminos de las matrices
% vect_aux2 = dv_e * (dV_dr' * despl);
% Dif11_Ax = vect_aux2(1);
% Dif11_Ay = vect_aux2(2);
% Dif11_Az = vect_aux2(3);
% 
% 
% % Duodecimo termino
% vect_aux = dv_r_red(1) * dAx_dv + dv_r_red(2) * dAy_dv + dv_r_red(3) * dAz_dv;
% Dif12_V = dot( dv_e, vect_aux);
% 
% vect_aux2 = dv_e * (dV_dv' * dv_r_red);     % Valdra 0
% Dif12_Ax = vect_aux2(1);
% Dif12_Ay = vect_aux2(2);
% Dif12_Az = vect_aux2(3);


% Por fin, ahora juntamos todos los terminos
% Dif_V = Dif1_V + Dif2_V + Dif3_V + Dif4_V + Dif5_V + Dif6_V + Dif7_V + Dif8_V + Dif9_V + Dif10_V + Dif11_V + Dif12_V;
% Dif_Ax = Dif1_Ax + Dif2_Ax + Dif3_Ax + Dif4_Ax + Dif5_Ax + Dif6_Ax + Dif7_Ax + Dif8_Ax + Dif9_Ax + Dif10_Ax + Dif11_Ax + Dif12_Ax;
% Dif_Ay = Dif1_Ay + Dif2_Ay + Dif3_Ay + Dif4_Ay + Dif5_Ay + Dif6_Ay + Dif7_Ay + Dif8_Ay + Dif9_Ay + Dif10_Ay + Dif11_Ay + Dif12_Ay;
% Dif_Az = Dif1_Az + Dif2_Az + Dif3_Az + Dif4_Az + Dif5_Az + Dif6_Az + Dif7_Az + Dif8_Az + Dif9_Az + Dif10_Az + Dif11_Az + Dif12_Az;

Dif_V = Dif1_V + Dif2_V + Dif3_V + Dif4_V + Dif5_V + Dif6_V;
Dif_Ax = Dif1_Ax + Dif2_Ax + Dif3_Ax + Dif4_Ax + Dif5_Ax + Dif6_Ax;
Dif_Ay = Dif1_Ay + Dif2_Ay + Dif3_Ay + Dif4_Ay + Dif5_Ay + Dif6_Ay;
Dif_Az = Dif1_Az + Dif2_Az + Dif3_Az + Dif4_Az + Dif5_Az + Dif6_Az;



Deriv = [Dif_V, Dif_Ax, Dif_Ay, Dif_Az];