% No es que la "v1" no funcione, pero este Matlab intenta obtener el mismo
% resultado usando sólo rotaciones de Wigner sin necesitar al SRIC,
% pensando que si alineamos la partícula "A" con el observador "B", la
% aceleración evaluada en A debería estar rotada respecto al valor que
% procede del SRIC.

% 20.09. Este enfoque también ha funcionado! Y con la prueba de recuperar
% v_ca_calc ha derrotado al enfoque de V1.

clear;

R = 1;
% Esta distancia está vista desde `U`.
num_iters = 10000;

v = 0.4;
% Prueba de control:
%v_cb  = [0, 0, 0];
% Prueba de interés:
%v_cb  = [0, 0.5, 0];

% Prueba más dura (tambien va bien 05.09.2026):
v_cb  = [0.15, -0.1, 0.5];

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

% Tiene que estar alineada con `u` para que la transformación de
% velocidades inicial no lo joda todo! --> Parece que no hace falta
% 19.09. Había -v_cb ¡Me cago en mi vida!
%v_ba_ini = Vel_Addition_Law(v_ca_ini, -v_cb)
v_ba_ini = Vel_Addition_Law(v_ca_ini, v_cb)
v_ba = v_ba_ini;

iter = 0;
while iter < num_iters
    
    if iter == 2500
        xxx = norm(v_ba);
        v_ca_calc = Vel_Addition_Law(v_ba, -v_cb)
        norm_v_ca_calc = norm(v_ca_calc)
    end
    
    a = -a_abs * [cos(fase), sin(fase), 0];
    dv = a * dt_prim;  
    
    % Puedo sacar la rotación de Wigner con v_cb que nunca cambia y v_ca que
    % es de fiar, esto es visto desde "C", pero mientras las velocidades
    % sean correctas y no uses -v_girada, ya está bien.
    v_ca = v * [-sin(fase), cos(fase), 0];
    RW_4 = RotacionWigner(v_ca, v_cb);
    RW_3 = RW_4(2:4, 2:4);
    
    % Aceleración vista desde el SRI del observador "B"
    dv_rot = dv * RW_3;
    
    % Ahora hacer un boost a la velocidad con la que `A` ve a `U` pero
    % aplicando la rotación
    x_mini_boost = Boost([1, 0, 0, 0], -dv_rot);  
    
    v_au = -v_ba;
    
    % El boost no necesita rotación ahora
    x_vel = Boost(x_mini_boost, v_au);
    
    v_ba_new = x_vel(2:4) / x_vel(1);    
    
    v_ba = v_ba_new;
    
    fase = fase + dif_fase;
    iter = iter + 1; 
end

fase_final = fase - 2 * pi;
%precesion_total = precesion_acum
v_ba_final  = v_ba
