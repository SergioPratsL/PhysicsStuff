% Hasta este punto, este script calcula la velocidad de la misma forma que
% el V1.

clear;

% addpath('C:\Fisica\Scripts matlab\HacerTuTeoria\8_origenes')

R = 1;
% Esta distancia está vista desde `U`.
num_iters = 40000;

v = 0.4;
% Prueba de control:
%v_cb  = [0, 0, 0];
% Prueba de interés:
%v_cb  = [0, 0.5, 0];

% Prueba más dura (tambien va bien):
v_cb  = [0.25, -0.2, 0.4];

gamma = fGamma(v);

% Frecuecia angular en el SRI local de la partícula en MCU
w_prim = gamma * v / R;

T_prim = 2 * pi / w_prim;
dt_prim = T_prim / num_iters;

fase_ini = 0;
fase = fase_ini;
dif_fase = 2 * pi / num_iters;

% w_prim^2 o gamma^2 * w^2 pero no ambas cosas!!
a_abs = w_prim^2 * R;

v_ca_ini = v * [-sin(fase_ini), cos(fase_ini), 0];
v_ac_ini = - v_ca_ini;

v_ab_ini = Vel_Addition_Law(v_cb, v_ca_ini)
v_ba_ini = Vel_Addition_Law(v_ca_ini, v_cb)
% Alineamos v_ab con B, así son C y A los que tienen que estar rotados
%v_ab = v_ab_ini;
v_ab = -v_ba_ini;

iter = 0;
precesion_acum = 0;

while iter < num_iters           
    if iter == num_iters / 4
        
        v_ca = v * [-sin(fase), cos(fase), 0]
        % RATIFICADO!! (Pero esto es el enfoque de v1).
        %v_cb_calc = Vel_Addition_Law(v_ab, -v_ca);
        %norm_v_cb_calc = norm(v_cb_calc);
        
        % -v_ab en verdad es v_ba, en el nuevo enfoque B mueve los hilos
        v_ca_calc = Vel_Addition_Law(-v_ab, -v_cb)
    end
    
    a = -a_abs * [cos(fase), sin(fase), 0];
    %dv = a * dt_prim;
    
    % Camino alternativo:
    % Enviar la velocidad a través de C --> B --> A sin rotar para que
    % llegue como siendo de B,  de esta forma es como si B está alineado
    % con C y con A, dejando la unión entre C y A como la que rota.
%     dp_c = a * dt_prim / gamma;
%     dp_c_4 = [0, dp_c];
%     dp_b_4 = Boost(dp_c_4, v_cb);
%     dp_4_alt = Boost(dp_b_4, -v_ab);       
    
    % Alternativa de la alternativa... o copia de la copia de la imitación
    E_c = a / gamma;
    B_c = [0,0,0];
    [E_b, B_b] = Boost_EM(E_c, B_c, v_cb);
    F_b = FLorentz(-v_ab, E_b, B_b);
    pot = dot(-v_ab, F_b);
    dp_b_4 = [pot, F_b] * dt_prim * fGamma(v_ab);
    dp_4 = Boost(dp_b_4, -v_ab);
    
    % Es dv porque en el SRI de A dp y dv son equivalentes (m=1)    
    dv = dp_4(2:4);    
  
    % Ahora hacer un boost a la velocidad con la que `A` ve a `U`.
    v_ab_new = Vel_Addition_Law(v_ab, dv);

    v_ca = v * [-sin(fase), cos(fase), 0];
    % Formula correcta cuando la aceleración viene directa del SRIC
    %RW = RotacionWigner(-v_ca, dv);
    % Quiero verificar que si saco la precesión que correspondería a v_ab,
    % el resultado falla.
    RW = RotacionWigner(v_ab, dv);
    
    RW_3 = RW(2:4, 2:4);
    
    % Esto va cuando la aceleración viene directa del SRIC
    v_ab_new = v_ab_new * RW_3';
    
    v_ab = v_ab_new;
    
    fase = fase + dif_fase;
    iter = iter + 1; 
end

fase_final = fase - 2 * pi;
v_ab_final  = v_ab
