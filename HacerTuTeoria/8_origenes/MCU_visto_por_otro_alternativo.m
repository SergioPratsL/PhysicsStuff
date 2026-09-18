% La idea de este Matlab es complementar al MCU_visto_por_otro.m haciendo
% que aquí la variación de la velocidad no la calcula el objeto en
% aceleración y se la reporta al "otro" sino que la calcula el "otro"
% directamente, con lo que las rotaciones a realizar serán diferentes.

clear;

R = 1;
% Esta distancia está vista desde `U`.
num_iters = 10000;

v = 0.4;
% Prueba de control:
%v_u  = [0, 0, 0];
% Prueba de interés:
%v_u  = [0, 0.5, 0];

% Prueba máss dura (tambien va bien 05.09.2026):
v_u  = [0.15, 0, 0.5];

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
% velocidades inicial no lo joda todo!
v_ua_ini = Vel_Addition_Law(v_ca_ini, -v_u)
v_ua = v_ua_ini;

iter = 0;
precesion_acum = 0;

% Obtenido de la prueba de control, para debug (obviamnete si v cambia,
% este valor ha de cambiar.
precesion_total_centro = 0.572332;
angulo_chungo = precesion_total_centro / num_iters;

% Esto se acerca, pero no llega :(
% while iter < num_iters
%     
%     if iter == 5000
%         abd = 2131;
%     end
%     
%     a = -a_abs * [cos(fase), sin(fase), 0];
%     dv = a * dt_prim;
%     
%     [dir_rot, angulo_W_T] = Angulo_Wigner(-v_ua, dv);
%         
%     % Ahora hacer un boost a la velocidad con la que `A` ve a `U` pero
%     % aplicando la rotación
%     x_mini_boost = Boost([1, 0, 0, 0], -dv); 
%     x_vel = Boost(x_mini_boost, -v_ua);    
%     v_ua_new = x_vel(2:4) / x_vel(1);        
%     
%     angulo = (angulo_W_T - angulo_chungo);
%     mat_rot = MatrizRotacion(dir_rot, angulo);
%     
%     precesion_acum = precesion_acum + angulo;
%     
%     v_ua_new = v_ua_new * mat_rot;
%     
%     v_ua = v_ua_new;
%     
%     fase = fase + dif_fase;
%     iter = iter + 1; 
% end


eje_chungo = [0, 0, 1];
mat_rot_1 = MatrizRotacion(eje_chungo, angulo_chungo); 
% Invirtiendo angulo_chungo y angulo_W_T no va...
%mat_rot_1 = MatrizRotacion(eje_chungo, -angulo_chungo); 

while iter < num_iters
    
    if iter == 2500
        abd = 2131;
    end
    
    a = -a_abs * [cos(fase), sin(fase), 0];
    dv = a * dt_prim;
        
    % Ahora hacer un boost a la velocidad con la que `A` ve a `U` pero
    % aplicando la rotación
    x_mini_boost = Boost([1, 0, 0, 0], -dv); 
    v_au_rot = -v_ua * mat_rot_1;
    
    x_vel = Boost(x_mini_boost, v_au_rot);    
    v_ua_new = x_vel(2:4) / x_vel(1);        
    
    [dir_rot, angulo_W_T] = Angulo_Wigner(-v_ua, dv);
    mat_rot_2 = MatrizRotacion(dir_rot, angulo_W_T);
    % Invirtiendo angulo_chungo y angulo_W_T no va...
    %mat_rot_2 = MatrizRotacion(dir_rot, -angulo_W_T);
    
    v_ua_new = v_ua_new * mat_rot_2;
    
    precesion_acum = precesion_acum + (angulo_W_T - angulo_chungo);
    
    v_ua = v_ua_new;
    
    fase = fase + dif_fase;
    iter = iter + 1; 
end

fase_final = fase - 2 * pi
precesion_total = precesion_acum
v_ua_final  = v_ua
