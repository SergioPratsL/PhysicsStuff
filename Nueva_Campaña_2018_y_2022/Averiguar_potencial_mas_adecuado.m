% No hay muchas constantes en este matlab, la permitividad y permeabilidad
% están fuera, de forma que la fuerza es del palo q1*q2/r^2
% Velocidad de la luz normalizada.
% En este caso, para simplificar no se acelera, sólo quiero comparar
% diferentes potenciales con la energía transferida al campo.

clear

% Masas de las particulas, si son electrones lo asimilo a la unidad
%m1 = 1;
%m2 = 1;

% Carga de las partículas, debo asumir que sean iguales o en todo caso de
% signo negativo, deben ser menores que uno si quiero que la distancia de
% separación minima sea 1 sin obtener valores extremos
q1 = 1;     % Si masa = 1 convenia elegir valor mas pequeño como 0.04
q2 = 1;


pos_1_ini = [-500, -0.5, 0];
pos_2_ini = [500, 0.5, 0];

% Velocidades, no tienen que ser iguales pero siempre son en la X por
% convenio. Expresadas en relación a la velocidad de la luz
v1_ini = [0.02, 0, 0];
v2_ini = [-0.4, 0, 0];


%consid_rad = 1;


prod_q = q1 * q2;
q1_cuad = q1^2;
q2_cuad = q2^2;

% Para saber cuánto se mueven por unidad de tiempo, en un primer momento
dif_vel = norm(v1_ini - v2_ini );


d_ini = norm(pos_1_ini - pos_2_ini);

% Iteraciones por octava (es decir, hasta que se divida o multiplique
% por dos la distancia inicial
num_it_octava = 5000;                % Maximo usado: 24000 (Pero la variación en la Ec crece)

dr_iteracion = d_ini / num_it_octava / 2;

dt_iteracion = dr_iteracion / dif_vel;

% Estimar el tamaño de los vectores para no tener que recortarlos (memoria
% no debería faltar para un algoritmo como este)
num_octavas_estimado = 25;

tam_vector = num_octavas_estimado * num_it_octava * 2;


pos_1 = zeros(tam_vector, 3);
pos_2 = zeros(tam_vector, 3);
vel_1 = zeros(tam_vector, 3);
vel_2 = zeros(tam_vector, 3);
% Obtener también el momento a partir del cual sacaremos la velocidad
mom_1 = zeros(tam_vector, 3);
mom_2 = zeros(tam_vector, 3);
t_iter = zeros(tam_vector, 3);

% if consid_rad == 1
%     acel_1 = zeros(tam_vector, 3);
%     acel_2 = zeros(tam_vector, 3);
%     Var_Ener_Rad_1 = 0;
%     Var_Ener_Rad_2 = 0;
%     Ener_Larmor_1 = 0;
%     Ener_Larmor_2 = 0;
% end

pos_1(1,:) = pos_1_ini;
pos_2(1,:) = pos_2_ini;
vel_1(1,:) = v1_ini;
vel_2(1,:) = v2_ini;
% mom_1(1,:) = m1 / sqrt( 1 - norm(vel_1(1,:))^2) * vel_1(1,:);
% mom_2(1,:) = m2 / sqrt( 1 - norm(vel_2(1,:))^2) * vel_2(1,:);


iter = 2;       % Se incrementa al principio
iter_old = 1;   % Se incrementa
t = 0;
% Simplemente es para forzar que se calcule la posición retardada a partir
% de asumir velocidad constante.
tr1 = -1;
tr2 = -1;

dist = d_ini;
d_ini_original = d_ini;      % d_ini se actualiza en cada iteración!

d_min = 10; % Por debajo de este valor se sale del bucle.

Dif_Ec_1 = 0;
Dif_Ec_2 = 0;

V_ini_1 = 0;

while iter < tam_vector && dist < (d_ini_original + 1) && dist > d_min
    t = t + dt_iteracion;
    
    pos_1(iter,:) = pos_1(iter_old,:) + vel_1(iter_old,:) * dt_iteracion;
    pos_2(iter,:) = pos_2(iter_old,:) + vel_2(iter_old,:) * dt_iteracion;
  
    
% Obtener la distancia retardada despreciando Y porque en estas distancias
% el efecto es mínimo, este if tiene un impacto muy pequeño (ha de tenerlo)
    if tr1 < 0
        tr1 = t - dist / (1- norm(vel_2(1,:)));
        pos_ret1 = pos_2(iter,:) - vel_2(iter_old,:) * (t - tr1);
        
        if tr1 >= 0
            indx_pos_ret1 = 1;
        end
        vel_ret1 = vel_2(1,:);
        
% Algoritmo para determinar la posición retardada        
    else
       seguir = 1;
       while seguir == 1
          indx_pos_ret1_next = indx_pos_ret1 + 1;
          dist_aux1 = norm( pos_1(iter,:) - pos_2(indx_pos_ret1_next,:)) ;
          
% Vemos si la luz de la siguiente posición es posterior a la luz retardada
% o no, si lo es estamos en el segmento adecuado, si no, hay que avanzar.
          if ( t - t_iter(indx_pos_ret1_next) > dist_aux1 )
            indx_pos_ret1 = indx_pos_ret1 + 1;  
          else
            seguir = 0;  
          end        
       end
       interval = t_iter(indx_pos_ret1_next) - t_iter(indx_pos_ret1);
       % cuanto más pequeño es cof, más cerca estamos de pos_ret1, por ello
       % debe ir multiplicado por (1-cof), no por cof!
       cof = abs(t - t_iter(indx_pos_ret1) - dist_aux1) / interval;
       
       % Hacemos interpolación lineal de la posición origen y de la velocidad
       pos_ret1 = (1-cof) * pos_2(indx_pos_ret1,:) + cof * pos_2(indx_pos_ret1_next,:); 
       vel_ret1 = (1-cof) * vel_2(indx_pos_ret1,:) + cof * vel_2(indx_pos_ret1_next,:);
       
    end
    
    
    if tr2 < 0
        tr2 = t - dist / (1- norm(vel_1(1,:)));
        pos_ret2 = pos_1(iter,:) - vel_1(iter_old,:) * (t - tr2);
        
        if tr2 >= 0
            indx_pos_ret2 = 1;
        end
        vel_ret2 = vel_1(1,:);
        
        %if consid_rad == 1
        %    acel_ret_2 = [0,0,0];  % Despreciar la aceleración retardada si lejos
        %end        
        
% Algoritmo para determinar la posición retardada        
    else
       seguir = 1;
       while seguir == 1
          indx_pos_ret2_next = indx_pos_ret2 + 1;
          dist_aux2 = norm( pos_2(iter,:) - pos_1(indx_pos_ret2_next,:)) ;
          
% Vemos si la luz de la siguiente posición es posterior a la luz retardada
% o no, si lo es estamos en el segmento adecuado, si no, hay que avanzar.
          if ( t - t_iter(indx_pos_ret2_next) > dist_aux2 )
            indx_pos_ret2 = indx_pos_ret2 + 1;  
          else
            seguir = 0;  
          end        
       end
       interval = t_iter(indx_pos_ret2_next) - t_iter(indx_pos_ret2);
       cof = abs(t - t_iter(indx_pos_ret2) - dist_aux2) / interval;
       
       % Hacemos interpolación lineal de la posición origen y de la velocidad
       pos_ret2 = (1-cof) * pos_1(indx_pos_ret2,:) + cof * pos_1(indx_pos_ret2_next,:);
       vel_ret2 = (1-cof) * vel_1(indx_pos_ret2,:) + cof * vel_1(indx_pos_ret2_next,:);

    end    
    
% Sacar el campo EM y la fuerza de Lorentz
    vect_retard_1 = pos_1(iter,:) - pos_ret1;
    vect_retard_2 = pos_2(iter,:) - pos_ret2;
    
    [E1, B1] = CampoInducido_sin_unidades(vect_retard_1, vel_ret1);
    [E2, B2] = CampoInducido_sin_unidades(vect_retard_2, vel_ret2);    
    
    F1 = FLorentz( vel_1(iter,:), E1, B1);
    F2 = FLorentz( vel_2(iter,:), E2, B2);   
    
    vel_1(iter,:) = v1_ini;
    vel_2(iter,:) = v2_ini;
    
    Dif_Ec_1 = Dif_Ec_1 + dot(F1, vel_1(iter,:)) * dt_iteracion;
    Dif_Ec_2 = Dif_Ec_2 + dot(F2, vel_2(iter,:)) * dt_iteracion;
    
    
    [V_1, A_1] = Potencial_LW( (pos_1(iter,:)-pos_ret1), vel_2(iter,:));
    V_1 = V_1 * prod_q;
    A_1 = A_1 * prod_q;
    Lagr_1 = V_1 + dot(A_1, vel_1(iter,:));     % signo deberia ser menos pero con mas va bien...

    [V_2, A_2] = Potencial_LW( (pos_2(iter,:)-pos_ret2), vel_1(iter,:));
    V_2 = V_2 * prod_q;
    A_2 = A_2 * prod_q;
    Lagr_2 = V_2 + dot(A_2, vel_2(iter,:));    
    
    if V_ini_1 == 0
        V_ini_1 = V_1;
        Lagr_ini_1 = Lagr_1;
        V_ini_2 = V_2;
        Lagr_ini_2 = Lagr_2;
    end    
    
    
% Por último evaluar si hay que cambiar de octava.
    d_act = norm(pos_1(iter,:) - pos_2(iter,:));
    
    if d_act <= (d_ini / 2) || d_act >= (d_ini * 2)
        d_ini = d_act;
        dr_iteracion = d_ini / num_it_octava / 2;
        dif_vel = norm(vel_1(iter,:) - vel_2(iter,:));
        dt_iteracion = dr_iteracion / dif_vel;  
    end
    
    
    iter_old = iter;
    t_iter(iter) = t; 
    dist = norm( pos_1(iter,:) - pos_2(iter,:));
    iter = iter + 1;
    
end


% Mostrar valores finales:
numero_de_iteraciones = iter_old

pos_1_final = pos_1(iter_old,:)
pos_2_final = pos_2(iter_old,:)

Dif_pot_1 = V_1 - V_ini_1
Dif_Lagr_1 = Lagr_1 - Lagr_ini_1
Dif_Ec_1 = Dif_Ec_1

Dif_pot_2 = V_2 - V_ini_2
Dif_Lagr_2 = Lagr_2 - Lagr_ini_2
Dif_Ec_2 = Dif_Ec_2

Dif_Ec_tot = Dif_Ec_1 + Dif_Ec_2


