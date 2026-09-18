% Este Matlab es la continuación de MCU_visto_por_otro.m y
% MCU_visto_por_otro_alternativo.m y su objetivo es ver cómo ha de cambiar
% la forma de transformar el tensor EM para que la fuerza de Lorentz lleve
% al nuevo valor de dp que se obtiene al cambiar la aceleración que ve el
% observador que no es el centro, el llamado "B".

clear;

R = 1;
fact_amp = 1;
num_iters = 10000 * fact_amp;

v = 0.4;
% Prueba de control:
%v_cb  = [0, 0, 0];

% Prueba de interés:
v_cb  = [0, 0.5, 0];

% Prueba algo más dura:
%v_cb  = [0, 0, 0.5];


% Prueba más dura (tambien va bien 10.09.2026):
%v_cb  = [0.15, 0, 0.5];

gamma = fGamma(v);

w_sric = v / R;
% Frecuecia angular en el SRI local de la partícula en MCU
w_prim = gamma * w_sric;

T_prim = 2 * pi / w_prim;
dt_prim = T_prim / num_iters;

fase_ini = 0;
fase = fase_ini;
dif_fase = 2 * pi / num_iters;

% w_prim^2 o gamma^2 * w^2 pero no ambas cosas!!
a_abs = w_prim^2 * R;
a_sric = a_abs / gamma;

v_ca_ini = v * [-sin(fase_ini), cos(fase_ini), 0];

% Tiene que estar alineada con `u` para que la transformación de
% velocidades inicial no lo joda todo!
v_ba_ini = Vel_Addition_Law(v_ca_ini, -v_cb)
v_ba = v_ba_ini;

iter = 0;
precesion_acum = 0;

% Este enfoque es un pocu tonto, me paro en un momento del proceso y
% experimento los valores que hay, el resto lo descarto. Ineficar, pero
% parece que de capacidad computacional voy sobrado.
iter_control = 3600 * fact_amp;

% Obtenido de la prueba de control, para debug (obviamnete si v cambia,
% este valor ha de cambiar.
precesion_total_centro = 0.572332;
angulo_sric = precesion_total_centro / num_iters;
eje_sric = [0, 0, 1];
mat_rot_1 = MatrizRotacion(eje_sric, angulo_sric); 

while iter < num_iters   
    a = -a_abs * [cos(fase), sin(fase), 0];
    dv = a * dt_prim;
        
    % Ahora hacer un boost a la velocidad con la que `A` ve a `U` pero
    % aplicando la rotación
    x_mini_boost = Boost([1, 0, 0, 0], -dv); 
    v_ab_rot = -v_ba * mat_rot_1;
    
    x_vel = Boost(x_mini_boost, v_ab_rot);    
    v_ba_new = x_vel(2:4) / x_vel(1);        
    
    [dir_rot, angulo_W_T] = Angulo_Wigner(-v_ba, dv);
    mat_rot_2 = MatrizRotacion(dir_rot, angulo_W_T);
    
    v_ba_new = v_ba_new * mat_rot_2;
    
    precesion_acum = precesion_acum + (angulo_W_T - angulo_sric);
    
    % Aquí paramos y hacemos la verdadera prueba!
    % Lo que veo es que aplicar una precesión que debería ser infinitesimal
    % al campo EM hará que ese campo casi no cambio, la fuerza cambiará
    % infinitesimalmente y eso no vale nada, hay que hacer la rotación de
    % Wigner entre los 3 sujetos...
    
    % PERO PARECE QUE LAS PRECESIONES DE WIGNER ENTRE LOS TRES SISTEMAS NO
    % DIERON LUGAR A NADA BUENO, ESTA VIA SE CANCELA
    if iter == iter_control        
       p_new = fGamma(v_ba_new) * v_ba_new;
       p_old = fGamma(v_ba) * v_ba;
       dp = p_new - p_old
       
       vect_clasico = Boost(x_mini_boost, -v_ba);
       v_ba_new_clasico = vect_clasico(2:4) / vect_clasico(1);
       p_new_clasico = fGamma(v_ba_new_clasico) * v_ba_new_clasico;
       dp_clasico = p_new_clasico - p_old       
     
       v_ca = v * [-sin(fase), cos(fase), 0];
       
       gamma_ba = fGamma(v_ba);
       
       % 15.09. También debo ver cómo transforma el dv visto en A a B con
       % las rotaciones y las pollas varias...
       dp_a = [0, dv];
       % No hace falta rotar, lo infinitesimal no hará nada... da el
       %  resultado clasico
       dp_boost = Boost(dp_a, -v_ba)
       
       
       
       % Usar v_ba_new que es la que corresponde a esta iteracion
       %[dir_rot_EM, angulo_EM] = Angulo_Wigner(-v_ba_new, -v_ca);
       % Esto es una prueba, veo que da un valor diferente... pues quiero
       % que sea el sistema central el que parta la pana.
       [dir_rot_EM, angulo_EM] = Angulo_Wigner(v_cb, v_ca);
       mat_rot_EM_4_4 = MatrizRotacion_4_4(dir_rot_EM, angulo_EM);
       mat_rot_EM = MatrizRotacion(dir_rot_EM, angulo_EM);
       
       % v_something es -v_cb, tiene un pequeño error que con más
       % iteraciones parece caer proporcionalmente. 
       %v_something = Vel_Addition_Law(-v_ba, -v_ca)
       %v_something2 = Vel_Addition_Law(-v_ca, -v_ba)       
       
       % Obtendré el campo EM visto desde el SRIC en un momento dado.
       E_sric = - a_sric * [cos(fase), sin(fase), 0];
       B_sric = [0,0,0];
       
       dp_sri = [0, E_sric * dt_prim * gamma];
       %dp_sri_rot = dp_sri * mat_rot_EM_4_4;
       dp_sri_rot = dp_sri;
       %v_cb_rot = v_cb * mat_rot_EM;
       v_cb_rot = v_cb;
       
       % Aplicando las rotaciones al destino pero no al origen, se obtiene
       % el valor dp_clasico.
       dp_u_from_sri_rot = Boost(dp_sri_rot, v_cb_rot);
       dp_u_from_sri = dp_u_from_sri_rot * mat_rot_EM_4_4';
       dp_from_sri = dp_u_from_sri(2:4)
       
       
       % Otro enfoque, B coge el campo EM de C y luego lo aplica a A pero
       % debe haber rotaciones que suceden en el SRI de B.

        [E_B, B_B] = Boost_EM(E_sric, B_sric, v_cb);
        %E_B_rot = E_B * mat_rot_EM;  
        E_B_rot = E_B;
        %B_B_rot = B_B * mat_rot_EM;
        B_B_rot = B_B;
        [E_Av2, B_Av2] = Boost_EM(E_B_rot, B_B_rot, v_ba);
        E_Av2 = E_Av2 * mat_rot_EM
        B_Av2 = B_Av2 * mat_rot_EM
       
       
       % parece que no voy a pasar de aquí.
% %        v_ba_rot = v_ba_new * mat_rot_EM;
% %        %v_ba_rot = v_ba_new;
% %        fB_v3 = FLorentz(v_ba_rot, E_B_rot, B_B_rot);
% %        
% %        dp_em_v3 = fB_v3 * dt_prim * gamma_ba
% %        norm_dp_em_v3 = norm(dp_em_v3)
       
       
       % Asumo que desde el SRIC puedes transformar los campos con un
       % simple boost "de toda la vida"
       %[E_A, B_A] = Boost_EM(E_sric, B_sric, v_ca_rot);
       [E_A, B_A] = Boost_EM(E_sric, B_sric, v_ca)
       
       F_A = TensorEM(E_A, B_A);
       
       F_A_rot = RotacionTensor(F_A, dir_rot_EM, angulo_EM);
       %F_A_rot = F_A;
       
       F_B_rot = BoostTensorEE(F_A_rot, -v_ba_new);
       
       % Para asegurar que tanto vba como el campo tras el boost rotan de
       % la misma manera (estoy haciendo muchas pruebas).
       angulo_retorno = angulo_EM;
       F_B = RotacionTensor(F_B_rot, dir_rot_EM, angulo_retorno);
       %F_B = F_B_rot;
       
       [E_B, B_B] = ExtraeEM_deTensor(F_B);
       
       mat_rot_EM = MatrizRotacion(dir_rot_EM, -angulo_retorno);
       v_ba_rot = v_ba_new * mat_rot_EM;
       
       fB = FLorentz(v_ba_rot, E_B, B_B);
       
       dp_em = fB * dt_prim * gamma_ba
       norm_dp_em = norm(dp_em)
       
       F_sric = TensorEM(E_sric, B_sric);
       F_sric_rot = RotacionTensor(F_sric, dir_rot_EM, angulo_EM);
       
       mat_rot_ad_hoc = MatrizRotacion(dir_rot_EM, angulo_EM);
       v_cb_rot = v_cb * mat_rot_ad_hoc;
       
       F_B_altern = BoostTensorEE(F_sric_rot, -v_cb_rot);
       F_B_altern = RotacionTensor(F_B_altern, dir_rot_EM, -angulo_EM);
       
       [E_B_altern, B_B_altern] = ExtraeEM_deTensor(F_B_altern);
       
       %[E_B_altern, B_B_altern] = Boost_EM(E_sric, B_sric, v_cb);
       fB_altern = FLorentz(v_ba, E_B_altern, B_B_altern);
       dp_em_altern = fB_altern * dt_prim * gamma_ba;
       norm_dp_altern = norm(dp_em_altern);
    end    
    
    
    v_ba = v_ba_new;
    
    fase = fase + dif_fase;
    iter = iter + 1; 
end

fase_final = fase - 2 * pi;
precesion_total = precesion_acum;
v_ba_final  = v_ba;
