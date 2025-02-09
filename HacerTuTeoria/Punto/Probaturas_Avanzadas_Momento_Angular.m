% Este script pretendía enmendar el Probaturas_Momento_Angular hecho en 
% 20.04.2023, lo que hice era un poco sui-generis y me parece que lo 
% entendí mal, parece que el momento angular no se conserva en un MCU pero
% si formo el vector P-L espero que eso si pueda devolver el L original

% En todas las pruebase se asume que R es 1 y además la masa en reposo
% del objeto que orbita también es 1

clear;

% Prueba 1 el boost es perpendicular a la distancia y al momento angular,
% ser paralelo a la distancia implica alterar la simultaneidad.
% u = [0.4, 0, 0];
% Gamma_lab = fGamma(u);

% x1 = [0, 1, 0];
% v1 = [0.3, 0, 0];
% x2 = [0,-1, 0];
% v2 = [-0.3, 0, 0];
% 
% v1_lab = Vel_Addition_Law(v1, -u);
% v2_lab = Vel_Addition_Law(v2, -u);
% 
% % Se llama Gamma, pero la función es fGamma... estoy vendido :P
% Sigma1_CM = fGamma(v1);
% Sigma2_CM = fGamma(v2);
% Sigma1_lab = fGamma(v1_lab);
% Sigma2_lab = fGamma(v2_lab);
% 
% 
% % No hay que modificar x1 ni x2 por ser perpendiculares a u.
% L1_CM = Sigma1_CM * cross(x1, v1)
% L2_CM = Sigma2_CM * cross(x2, v2)
% L_CM__Prueba1 = L1_CM + L2_CM
% 
% L1_lab = Sigma1 * cross(x1, v1_lab)
% L2_lab = Sigma2 * cross(x2, v2_lab)
% L_lab_Prueba1 = L1_lab + L2_lab
% 
% ratio_lab_CM_Prueba1 = norm(L_lab_Prueba1) / norm(L_CM_Prueba1)
% Parece que el momento angular crece Sigma

% Prueba 2, ahora el boost es en Z, la dirección del momento angular
% u = [0, 0, 0.4];
% Gamma_lab = fGamma(u);
% 
% x1 = [-1, 0, 0];
% v1 = [0, 0.3, 0];
% x2 = [1, 0, 0];
% v2 = [0, -0.3, 0];
% 
% v1_lab = Vel_Addition_Law(v1, -u);
% v2_lab = Vel_Addition_Law(v2, -u);
% 
% Sigma1_CM = fGamma(v1);
% Sigma2_CM = fGamma(v2);
% Sigma1_lab = fGamma(v1_lab);
% Sigma2_lab = fGamma(v2_lab);
% 
% L1_CM = Sigma1_CM * cross(x1, v1)
% L2_CM = Sigma2_CM * cross(x2, v2)
% L_CM_Prueba2 = L1_CM + L2_CM
% 
% % Dividir por Gamma_lab es efecto de la contracción de longitudes.
% L1_lab = Sigma1_lab * cross(x1 / Gamma_lab, v1_lab)
% L2_lab = Sigma2_lab * cross(x2 / Gamma_lab, v2_lab)
% L_lab_Prueba2 = L1_lab + L2_lab
% 
% ratio_lab_CM_Prueba2 = norm(L_lab_Prueba2) / norm(L_CM_Prueba2)


% Prueba 3, lo mismo pero cambiando la posición que el boost sea
% paralelo a la posición de manera que hay destiempo
% El CM se mueva a u respecto al laboratorio.
% u = [0.4, 0, 0];
% Gamma_lab = fGamma(u);
% R = 1;
% v = 0.4;
% 
% x1 = [R, 0, 0];
% v1 = [0, v, 0];
% x2 = [-R, 0, 0];
% v2 = [0, -v, 0];
% 
% Sigma1_CM = fGamma(v1);
% Sigma2_CM = fGamma(v2);
% 
% L1_CM = Sigma1_CM * cross(x1, v1)
% L2_CM = Sigma2_CM * cross(x2, v2)
% L_CM_Prueba3 = L1_CM + L2_CM
% 
% % No hay simultaneidad entre CM y lab, por ello el CM verá las dos masa
% % haciendo un ángulo menor que 180 grados (y el complementario mayor que
% % 180º, obviamente, además, las distancias del "circuito" por el que va la
% % onda se contraen en la dirección de v...
% 
% % La diferencia de fase depende del destiempo, que depende de u y de la
% % velocidad angular, que depende de v.
% %dif_fase = -v * norm(u);
% dif_fase_1 = PosicionPresenteLab_MCU(0, v, u);
% dif_fase_2 = PosicionPresenteLab_MCU(pi, v, u);
% 
% x1_CM_destiempo = [cos(dif_fase_1), sin(dif_fase_1), 0];
% v1_CM_destiempo = v * [-sin(dif_fase_1), cos(dif_fase_1), 0];
% 
% % La diferencia de fase se invierte entre 1 y 2 pues su posición es opuesta
% x2_CM_destiempo = [cos(pi + dif_fase_2), sin(pi + dif_fase_2), 0];
% v2_CM_destiempo = v * [-sin(pi + dif_fase_2), cos(pi - dif_fase_2), 0];
% 
% x1_lab = ContraeDireccion(x1_CM_destiempo, u);
% x2_lab = ContraeDireccion(x2_CM_destiempo, u);
% 
% v1_lab = Vel_Addition_Law(v1_CM_destiempo, -u);
% v2_lab = Vel_Addition_Law(v2_CM_destiempo, -u);
% 
% Sigma1_lab = fGamma(v1_lab);
% Sigma2_lab = fGamma(v2_lab);
% 
% L1_lab = Sigma1_lab * cross(x1_lab, v1_lab)
% L2_lab = Sigma2_lab * cross(x2_lab, v2_lab)
% L_lab_Prueba3 = L1_lab + L2_lab
% 
% ratio_lab_CM_Prueba3 = norm(L_lab_Prueba3) / norm(L_CM_Prueba3)
% % ratio_lab_CM_Prueba3 =   1.0919 (con v=0.3).

% Esto no va...



% ------------------------    PARTE 2    ------------------------
% Voy a asumir un MCU con un sistema de dos masas, una "infinita"
% en el centro que no se mueve y otra que se mueve, toda contribución a L 
% proviene de la parte movil pero el CM es donde está la masa infinita
% al cambiar de SRI es obvio que el momento angular cambiará, no hace 
% falta considerar la "PosicionPresenteLab_MCU", sólo hay que asumir que
% ves que la partícula está en una posiciónm transformas velocidades,
% contraes longitudes y ya está...

% Prueba 4: varias fases pero u = 0.4x, por tanto estamos en plano XY.
% u = [0.4, 0, 0];
% Gamma_lab_CM = fGamma(u);
% 
% %fase = 0;
% %fase = pi/2;
% %fase = pi/6;
% fase = pi;
% v_mcu = 0.3;
% Gamma_CM = fGamma(v_mcu);
% 
% % Posición en sistema MCU que ve la partícula
% x = [cos(fase), sin(fase), 0];
% v = v_mcu * [-sin(fase), cos(fase), 0];
% 
% L_CM = Gamma_CM * cross(x, v)
% 
% v_lab = Vel_Addition_Law(v, -u);
% Gamma_lab = fGamma(v_lab);
% 
% x_lab = ContraeDireccion(x, u);
% 
% L_lab = Gamma_lab * cross(x_lab, v_lab)
% 
% % Usar el momento de masa, no sé si t debe ser 0 pero aqui espero que sí.
% N_lab = Gamma_lab * x_lab;
% 
% M_lab = TensorMomAng(N_lab, L_lab);
% 
% P_lab = Momento_4(-u);
% 
% W_lab = Pseudovector_P_L(M_lab, P_lab)
% 
% % u es positivo por ser -(-u), es más fácil esto que hacer N_CM, J_CM y
% % recalcular W... de hecho bastaría con J_CM pero bueno...
% W_CM = Boost(W_lab, u)

% fase = 0; % W_CM = [0         0         0    0.3145] --> W0 = 0, W123 = L ¡Eureka!
% fase = pi/2; % W_CM = [0         0         0    0.3145] 
% fases pi/6 o pi, mismo resultado :)


% Prueba 5: fase no cambia, pero U puede tomar cualquier dirección.
% u_norm = 0.4;
% 
% %u = [0, 0, 1];
% u = [0, -1, 1];
% %u = [2, -1, 1];
% u = NormalizaCustom(u, u_norm);
% 
% Gamma_lab_CM = fGamma(u);
% 
% fase = 0;
% v_mcu = 0.3;
% Gamma_CM = fGamma(v_mcu);
% 
% % Posición en sistema MCU que ve la partícula
% x = [cos(fase), sin(fase), 0];
% v = v_mcu * [-sin(fase), cos(fase), 0];
% 
% L_CM = Gamma_CM * cross(x, v)
% 
% v_lab = Vel_Addition_Law(v, -u);
% Gamma_lab = fGamma(v_lab);
% 
% x_lab = ContraeDireccion(x, u);
% 
% L_lab = Gamma_lab * cross(x_lab, v_lab)
% 
% % Usar el momento de masa, no sé si t debe ser 0 pero aqui espero que sí.
% N_lab = Gamma_lab * x_lab;
% 
% M_lab = TensorMomAng(N_lab, L_lab);
% 
% P_lab = Momento_4(-u);
% 
% W_lab = Pseudovector_P_L(M_lab, P_lab)
% 
% % u es positivo por ser -(-u), es más fácil esto que hacer N_CM, J_CM y
% % recalcular W... de hecho bastaría con J_CM pero bueno...
% W_CM = Boost(W_lab, u)
% u = [0, 0, 0.4]; -> W_CM = -0.0000         0         0    0.3145
% u = 0   -0.2828    0.2828 -> W_CM = -0.0000         0   -0.0000    0.3145
% u = 0.3266   -0.1633    0.1633 --> W_CM = -0.0000    0.0000   -0.0000    0.3145

% Así las cosas parecen conservarse :).


% Prueba 6: busqueda del potencial vector...
 u_norm = 0.4;
 %u_norm = 0;
 v_mcu = 0.3;
 
 %u = NormalizaCustom([1, 0, 0], u_norm);
 %W_lab_1 = Vector_PL_otro_SRI(u, 0, v_mcu)
 %W_lab_2 = Vector_PL_otro_SRI(u, pi, v_mcu)
 %W_lab_3 = Vector_PL_otro_SRI(u, pi/3, v_mcu)
 % Mismo valor en los tres casos: [0         0         0    0.3145]

% u = NormalizaCustom([0, 1, 1], u_norm);
%
% W_lab_1 = Vector_PL_otro_SRI(u, 0, v_mcu);
% W_lab_2 = Vector_PL_otro_SRI(u, pi/4, v_mcu);
% W_lab_3 = Vector_PL_otro_SRI(u, pi/2, v_mcu);
% % Mismo valor en los tres casos: [0.0971         0    0.0143    0.3288]
 
 u = NormalizaCustom([-2, 1, 1], u_norm);
 
 %W_lab_1 = Vector_PL_otro_SRI(u, 0, v_mcu)
 %W_lab_2 = Vector_PL_otro_SRI(u, pi/4, v_mcu)
% W_lab_3 = Vector_PL_otro_SRI(u, pi/2, v_mcu)
% Mismo valor en los tres casos: [0.0560   -0.0095    0.0048    0.3193]


% Prueba 7: confirmar que el momento angular no cambiará si nos alejamos
% del sistema en un sistema con una masa mayor mucho mayor que la masa en
% rotación, esta es la prueba FINAL.

% Cada llamada a esta función es una prueba.
%PruebaFinalOrigenNoCM([1,0,0], 0, 0)
%PruebaFinalOrigenNoCM([1,0,0], 25, 0)
% W_tot = [0   0    0    0.3145]
% PruebaFinalOrigenNoCM([1,0,0], 0, pi/2)
% PruebaFinalOrigenNoCM([1,0,0], 25, pi/2)
% W_tot = [0   0    0    0.3145]
%PruebaFinalOrigenNoCM([1,0,0], 0, 0.75 * pi)
%PruebaFinalOrigenNoCM([1,0,0], 25, 0.75 * pi)
% W_tot = [0   0    0    0.3145]

%PruebaFinalOrigenNoCM([-2,1,1], 0, 0);
%PruebaFinalOrigenNoCM([-2,1,1], 40, 0);
%PruebaFinalOrigenNoCM([-2,1,1], 10, pi/3);
PruebaFinalOrigenNoCM([-2,1,1], -30, pi);
% W_tot = [0.0560   -0.0095    0.0048    0.3193]

%0.0560   -0.0017    0.0363    0.3035


% Victoria decisiva, pero me ha extrañado el que centro se ha de quedar
% la diferencia entre el momento total "si no hubiera MCU" menos el momento
% actual del objeto que orbita... yo pensaba que debería tenerse en cuenta
% que lo que es simultaneo en el CM no lo es en Lab... pero tal vez es que
% también se debe considerar el momento del potencial vector asignándolo al
% centro... raro.


function W_lab = Vector_PL_otro_SRI(u, fase, v_mcu)

x = [cos(fase), sin(fase), 0];
v = v_mcu * [-sin(fase), cos(fase), 0];

v_lab = Vel_Addition_Law(v, -u);
Gamma_lab = fGamma(v_lab);

x_lab = ContraeDireccion(x, u);

L_lab = Gamma_lab * cross(x_lab, v_lab)
N_lab = Gamma_lab * x_lab;
M_lab = TensorMomAng(N_lab, L_lab);

P_lab = Momento_4(-u);

W_lab = Pseudovector_P_L(M_lab, P_lab);

end

% Obtiene la velocidad ya en el SRI del lab
function v_lab = ObtenVelocFase(fase, v_mcu, u)
    v_ori = v_mcu * [-sin(fase), cos(fase), 0];
    v_lab = Vel_Addition_Law(v_ori, -u);
end

function x_lab = ObtenVectorFase(fase, u, offSet)
    x_CM = [cos(fase), sin(fase), 0];
    x_lab = ContraeDireccion(x_CM, u) + offSet;
end

function W = ObtenVectorP_L(L, N, u)
    M = TensorMomAng(N, L);
    P = Momento_4(-u);
    W = Pseudovector_P_L(M, P);
end

% Esto mostrará los resultados, es la prueba última de este script :).
function PruebaFinalOrigenNoCM(u_not_norm, t, fase)

    u_norm = 0.4;
    v_mcu = 0.3;
    Gamma_CM = fGamma(u_norm);

    m = 1;
    M = 10^6;       % Centro super pesado

    u = NormalizaCustom(u_not_norm, u_norm);

    % spoiler, no se conserva (oscila periódicamente)
    MasaEnergia_CM = M + fGamma(v_mcu) * m;  % Purista.
    p_total = Gamma_CM * u * MasaEnergia_CM;

    offSet = u * t;
    % Algo así cambia el W total, pero no depende de la fase (ni de t).
    %offSet = offSet + [0, 1, 2];
    
    % sat es para abreviar satélite
    x_sat = ObtenVectorFase(fase, u, offSet);
    v_sat = ObtenVelocFase(fase, v_mcu, u);
    Gamma_sat = fGamma(v_sat);

    p_sat = Gamma_sat * v_sat * m;

    L_sat = cross(x_sat, p_sat);
    N_sat = Gamma_sat * m * x_sat;

    W_sat = ObtenVectorP_L(L_sat, N_sat, u)

    % Necesito dif_fase para que el centro tenga el momento que debería
    % compensar el momento de la partícula, pero en una posición un poco
    % diferente por las no simultaneidad!
    % El signo menos porque el "origen de tiempo" es el satelite, no el centro.
    dif_fase = - PosicionPresenteLab_MCU(0, v_mcu, u);

    v_sat_dif_fase = ObtenVelocFase(fase + dif_fase, v_mcu, u);
    p_sat_dif_fase = fGamma(v_sat_dif_fase) * m * v_sat_dif_fase;

    % Se asume que si todo fuera inercial se debe cumplir que p_M+p_m= p_total
    %p_centro = p_total - p_sat_dif_fase;
    % No debería pero parece que va, supera todas las pruebas 
    % mientras que lo desfasado, no
    p_centro = p_total - p_sat;            

    E_centro = sqrt(M^2 + norm(p_centro)^2);

    L_centro = cross(offSet, p_centro);
    % No usar Gamma_CM * M sino E_centro
    N_centro = E_centro * offSet;

    W_centro = ObtenVectorP_L(L_centro, N_centro, u)
    
    W_tot = W_sat + W_centro

    % Sólo para comprobar que esto no se mantiene ni pa atrás.
    % Se mantiene con la distancia pero no con la fase.
    L_tot = L_centro + L_sat
end


% Esto finalmente no resulto de utilidad:
% function dif_fase = PosicionPresenteLab_MCU(fase, v, u)
%     % Algoritmo iterativo para encontrar la posición en la que estaria la
%     % partícula en el "presente" del SRI laboratio si en el presente del Centro
%     % de Masa se encuentra en la posición dada por la fase.
% 
%     % Fase se usa para determinar la posición en el "presente" del MCU y a
%     % partir de ella se determina v, que de momento no podra elegir signo, se
%     % asume MCU sobre el plano XY.
%     % El Radio del MCU se asume que es 1
%     % v es la velocidad del MCU en valor absoluto, c=1
%     % u es la velocidad a la que el Lab ve moverse al CM, es un vector
% 
%     x_CM = posicion_fase(fase);
% 
%     % El Gamma(u) se va porque debemos medir en tiempo propio del CM
%     dif_t  = - dot(u, x_CM);
% 
%     dif_fase_ini = dif_t * v;
%     error_fase = abs(dif_fase_ini);
% 
%     iter = 0;
% 
%     while error_fase > 0.0000000001 && iter < 50
%         x_t_despl = posicion_fase(fase + dif_fase_ini);
% 
%         dif_t  = - dot(u, x_t_despl);
%         dif_fase = dif_t * v;
% 
%         error_fase = abs(dif_fase - dif_fase_ini);
% 
%         dif_fase_ini = dif_fase;
%         iter = iter + 1;
%     end
% 
%     % Esto da cero, lo que indica que la lógica es correcta.
%     %t_final_test = fGamma(u) * ( dif_t + dot(u, x_t_despl))
% 
% end
% 
% function x = posicion_fase(fase)
%     x = [cos(fase), sin(fase), 0];

