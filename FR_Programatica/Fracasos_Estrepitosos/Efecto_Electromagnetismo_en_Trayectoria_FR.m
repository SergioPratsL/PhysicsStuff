% No hay muchas constantes en este matlab, la permitividad y permeabilidad
% están fuera, de forma que la fuerza es del palo q1*q2/r^2
% Velocidad de la luz normalizada.
% en es

clear

% Masas de las particulas, si son electrones lo asimilo a la unidad
% m1 = 1;
% m2 = 1;

% Carga de las partículas, debo asumir que sean iguales o en todo caso de
% signo negativo, deben ser menores que uno si quiero que la distancia de
% separación minima sea 1 sin obtener valores extremos
%q1 = 0.04;
%q2 = 0.04;
%q1 = 0.1;
%q2 = 0.1;

%pos_1_ini = [-500, -0.5, 0];
%pos_2_ini = [500, 0.5, 0];

% Velocidades, no tienen que ser iguales pero siempre son en la X por
% convenio. Expresadas en relación a la velocidad de la luz
%v1_ini = [0.1, 0, 0];
%v2_ini = [-0.1, 0, 0];

% Prueba 2: 1 es un electrón, 2 es un protón super pesado y quieto
%m1 = 1;
%m2 = 1000000;

%q1 = 0.04;
%q2 = -0.04;

%pos_1_ini = [-1000, 1, 0];
%pos_2_ini = [0,0,0];

%v1_ini = [0.2, 0, 0];
%v2_ini = [0, 0, 0];


% Prueba 3: la aparente paradoja.
% Aparente paradoja resuelta! Ciertamente la energía final es mayor que la
% energía final en los electrones que se acercan, para compensar la energia
% negativa del campo radiado (tras aplicar la pertinente correccion a la
% energía), cualquier parecido con la realidad es pura coincidencia pero la
% teoria es correcta
m1 = 1;
m2 = 1;

q1 = 1;
q2 = 1;     % A muerte!

% Choque frontal, viniendo desde muy lejos --> 23/05/2022 se introduce una
% separación.
%pos_1_ini = [-500000000, 0.5, 0];
%pos_2_ini = [500000000, -0.5, 0];

pos_1_ini = [-50000000, 0, 0];
pos_2_ini = [50000000, 0, 0];

% DEbug
dist_old = norm(pos_1_ini - pos_2_ini);

v1_ini = [0.6, 0, 0];
v2_ini = [-0.6, 0, 0];

% 15/07/2018. Joder! Me olvide de considerar el efecto del campo radiado
% sobre la otra partícula
consid_rad = 1;

consid_freact = 0;

if consid_freact == 1 && consid_rad ~= 1
    fprintf('¡¡¡consider_rad no puede ser 0 mientras que consid_freact es 1, no hay reacción sin radiación!!!')
    return
end

prod_q = q1 * q2;
q1_cuad = q1^2;
q2_cuad = q2^2;

% Para saber cuánto se mueven por unidad de tiempo, en un primer momento
dif_vel = norm(v1_ini - v2_ini );

dr_iteracion_min = 10000000;
dt_iteracion_min = 10000000;

d_ini = norm(pos_1_ini - pos_2_ini);

% Iteraciones por octava (es decir, hasta que se divida o multiplique
% por dos la distancia inicial
num_it_octava = 2000;                % Maximo usado: 24000 (Pero la variación en la Ec crece)

dr_iteracion = d_ini / num_it_octava / 2;

dt_iteracion = dr_iteracion / dif_vel;

% Estimar el tamaño de los vectores para no tener que recortarlos (memoria
% no debería faltar para un algoritmo como este)
%num_octavas_estimado = 25;
num_octavas_estimado = round(log2(norm(pos_2_ini - pos_1_ini))) + 10;      


%tam_vector = 2^(num_octavas_estimado) * num_it_octava  * 2;
tam_vector = num_octavas_estimado * num_it_octava * 2;


pos_1 = zeros(tam_vector, 3);
pos_2 = zeros(tam_vector, 3);
vel_1 = zeros(tam_vector, 3);
vel_2 = zeros(tam_vector, 3);
% Obtener también el momento a partir del cual sacaremos la velocidad
mom_1 = zeros(tam_vector, 3);
mom_2 = zeros(tam_vector, 3);
t_iter = zeros(tam_vector, 3);

if consid_rad == 1
    acel_1 = zeros(tam_vector, 3);
    acel_2 = zeros(tam_vector, 3);
    Var_Ener_Rad_1 = 0;
    Var_Ener_Rad_2 = 0;
    Ener_Larmor_1 = 0;
    Ener_Larmor_2 = 0;
end

% DEBUG
if consid_freact == 1
    dif_p_fr = zeros(tam_vector, 3);    
    dif_p_norm = zeros(tam_vector, 3);    
    
    dif_p_fr_1_old  = [0,0,0];    
    
end

pos_1(1,:) = pos_1_ini;
pos_2(1,:) = pos_2_ini;
vel_1(1,:) = v1_ini;
vel_2(1,:) = v2_ini;
mom_1(1,:) = m1 / sqrt( 1 - norm(vel_1(1,:))^2) * vel_1(1,:);
mom_2(1,:) = m2 / sqrt( 1 - norm(vel_2(1,:))^2) * vel_2(1,:);


iter = 2;       % Se incrementa al principio
iter_old = 1;   % Se incrementa
t = 0;
% Simplemente es para forzar que se calcule la posición retardada a partir
% de asumir velocidad constante.
tr1 = -1;
tr2 = -1;

dist = d_ini;
d_ini_original = d_ini;      % d_ini se actualiza en cada iteración!

v_max_1 = 0;
v_max_2 = 0;
v_min_1 = 1;
v_min_2 = 1;
v_abs_1 = norm(vel_1(1,:));
v_abs_2 = norm(vel_2(1,:));

V_ini_1 = 0;

while iter < tam_vector && dist < (d_ini_original + 1)
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
        
        if consid_rad == 1
            acel_ret_1 = [0,0,0];  % Despreciar la aceleración retardada si lejos
        end
     
% Potencial sólo se calcula sobre uno de los dos porque sólo pienso fijarme
% en el caso en que se acercan a la misma velocidad, es para entender de
% donde proceden las diferencias de energía
        if V_ini_1 == 0
            [V_ini_1, A_ini_1] = Potencial_LW( (pos_1(iter,:)-pos_ret1), v2_ini);
            V_ini_1 = V_ini_1 * prod_q;
            A_ini_1 = A_ini_1 * prod_q;
            Lagr_ini_1 = V_ini_1 - dot(A_ini_1, v1_ini);
        end
        
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
       
       if cof > 1
           adas = 1234; % Esto no ocurre ni en el caso desviado.
       end
       
       % Hacemos interpolación lineal de la posición origen y de la velocidad
       pos_ret1 = (1-cof) * pos_2(indx_pos_ret1,:) + cof * pos_2(indx_pos_ret1_next,:); 
       vel_ret1 = (1-cof) * vel_2(indx_pos_ret1,:) + cof * vel_2(indx_pos_ret1_next,:);
       if consid_rad == 1
           acel_ret_1 = (1-cof) * acel_2(indx_pos_ret1,:) + cof * acel_2(indx_pos_ret1_next,:);
       end
       
    end
    
    
    if tr2 < 0
        tr2 = t - dist / (1- norm(vel_1(1,:)));
        pos_ret2 = pos_1(iter,:) - vel_1(iter_old,:) * (t - tr2);
        
        if tr2 >= 0
            indx_pos_ret2 = 1;
        end
        vel_ret2 = vel_1(1,:);
        
        if consid_rad == 1
            acel_ret_2 = [0,0,0];  % Despreciar la aceleración retardada si lejos
        end        
        
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
       if consid_rad == 1
           acel_ret_2 = (1-cof) * acel_1(indx_pos_ret2,:) + cof * acel_1(indx_pos_ret2_next,:);
       end
       
    end    
    
% Sacar el campo EM y la fuerza de Lorentz
    vect_retard_1 = pos_1(iter,:) - pos_ret1;
    vect_retard_2 = pos_2(iter,:) - pos_ret2;
    
    [E1, B1] = CampoInducido_sin_unidades(vect_retard_1, vel_ret1);
    [E2, B2] = CampoInducido_sin_unidades(vect_retard_2, vel_ret2);    
    
    if consid_rad == 1
        [E_rad_1, B_rad_1] = CampoRadiadoCargaAcelerada2( vect_retard_1, vel_ret1, acel_ret_1);
        [E_rad_2, B_rad_2] = CampoRadiadoCargaAcelerada2( vect_retard_2, vel_ret2, acel_ret_2);      
        
        E1 = E1 + E_rad_1;
        B1 = B1 + B_rad_1;
        E2 = E2 + E_rad_2;
        B2 = B2 + B_rad_2;  
        
        % Calculo adicional a realizar: efecto del campo radiado del otro.
        %Var_Ener_Rad_1 = Var_Ener_Rad_1 + dot(E_rad_1, vel_1(iter_old, :)) * dt_iteracion * prod_q;
        %Var_Ener_Rad_2 = Var_Ener_Rad_2 + dot(E_rad_2, vel_2(iter_old, :)) * dt_iteracion * prod_q;
    end
    
    F1 = FLorentz( vel_1(iter,:), E1, B1);
    F2 = FLorentz( vel_2(iter,:), E2, B2);       
    
    dif_p_1 = F1 * dt_iteracion * prod_q;
    dif_p_2 = F2 * dt_iteracion * prod_q;
    
    % 24.05.2022. ¡Por fin la fuerza de reacción entra en juego!
    if consid_freact == 1 && iter > 5
        % dif_p_fr es igual a la fuerza de reacción que es proporcional a
        % 1/dt_interaccion multiplicado por el tiempo qu dura...
        dif_p_fr_1 = (2/3)*q1^2/m1 * (acel_1(iter-1,:) - acel_1(iter-2,:));
        dif_p_fr_2 = (2/3)*q2^2/m2 * (acel_2(iter-1,:) - acel_2(iter-2,:));
        
        % DEBUG
        dif_p_norm(iter, :) = dif_p_1;
        dif_p_fr(iter, :) = dif_p_fr_1;
        
        % DEBUG
        if abs(dif_p_fr_1) / abs(dif_p_1) > 1
            iter = iter;
        end
        
        % DEBUG
        xxx = sign(dif_p_fr_1(1)) * sign(dif_p_fr_1_old(1));
        if xxx == -1
            xxx = xxx;
        end      
        dif_p_fr_1_old = dif_p_fr_1;
            
        dif_p_1 = dif_p_1 + dif_p_fr_1;
        dif_p_2 = dif_p_2 + dif_p_fr_2;
    end
    
    mom_1(iter,:) = mom_1(iter_old,:) + dif_p_1;
    mom_2(iter,:) = mom_2(iter_old,:) + dif_p_2;
        
% v^2(1+p^2) = p^2 / m^2   
    mom_1_cuadr = mom_1(iter,1)^2 + mom_1(iter,2)^2 + mom_1(iter,3)^2;
    mom_2_cuadr = mom_2(iter,1)^2 + mom_2(iter,2)^2 + mom_2(iter,3)^2;
    v_abs_1 = sqrt( mom_1_cuadr / (1 + mom_1_cuadr) ) / m1;
    v_abs_2 = sqrt( mom_2_cuadr / (1 + mom_2_cuadr) ) / m2;
    
    vel_1(iter,:) = v_abs_1 * mom_1(iter,:) / norm(mom_1(iter,:));
    vel_2(iter,:) = v_abs_2 * mom_2(iter,:) / norm(mom_2(iter,:)); 
    % DEBUG
    %vel_1(iter,:) = v_abs_1 * [sign(mom_1(iter,1)),0,0];
    %vel_2(iter,:) = v_abs_2 * [sign(mom_2(iter,2)),0,0];
    
    if v_abs_1 > v_max_1 
        v_max_1 = v_abs_1;
    elseif v_abs_1 < v_min_1 
        v_min_1 = v_abs_1;
    end
    
    if v_abs_2 > v_max_2 
        v_max_2 = v_abs_2;
    elseif v_abs_2 < v_min_2 
        v_min_2 = v_abs_2;
    end    
    
    if consid_rad == 1 && iter > 2
        
        Sigma_1 = 1 / sqrt( 1 - v_abs_1^2);
        Sigma_2 = 1 / sqrt( 1 - v_abs_2^2);
        
        % Aproximación guarra de la aceleración medida en el laboratorio
        % servirá porque la velocidad es mayormente en componente X
        %acel_1(iter,:) = [dif_p_1(1) / Sigma_1, dif_p_1(2), dif_p_1(3)] / (m1 * Sigma_1 * dt_iteracion);
        %acel_2(iter,:) = [dif_p_2(1) / Sigma_2, dif_p_2(2), dif_p_2(3)] / (m2 * Sigma_2 * dt_iteracion);        
        
        acel_1(iter,:) = ( vel_1(iter,:) - vel_1((iter-1),:) ) / dt_iteracion;
        acel_2(iter,:) = ( vel_2(iter,:) - vel_2((iter-1),:) ) / dt_iteracion;         
        
        %Ener_Larmor_1 = Ener_Larmor_1 + norm(acel_1(iter,:))^2 * (2/3) * q1_cuad * dt_iteracion;
        %Ener_Larmor_2 = Ener_Larmor_2 + norm(acel_1(iter,:))^2 * (2/3) * q2_cuad * dt_iteracion;
        factor1 = (2/3) * Sigma_1^6 * (norm(acel_1(iter,:))^2 - norm(cross(vel_1(iter,:), acel_1(iter,:)))^2) * q1_cuad;
        factor2 = (2/3) * Sigma_2^6 * (norm(acel_2(iter,:))^2 - norm(cross(vel_2(iter,:), acel_2(iter,:)))^2) * q2_cuad;
        Ener_Larmor_1 = Ener_Larmor_1 + factor1 * dt_iteracion;
        Ener_Larmor_2 = Ener_Larmor_2 + factor2 * dt_iteracion;
    end    
    
    
% Por último evaluar si hay que cambiar de octava.
    dist = norm( pos_1(iter,:) - pos_2(iter,:));
    
    % DEbug
    if dist > dist_old
       dist = dist; 
    end
    dist_old = dist;
    
    if dist <= (d_ini / 2) || dist >= (d_ini * 2)
        d_ini = dist;
        dr_iteracion = d_ini / num_it_octava / 2;
        dif_vel = norm(vel_1(iter,:) - vel_2(iter,:));
        dt_iteracion = dr_iteracion / dif_vel
        
        %DEbug
        if( dr_iteracion < dr_iteracion_min)
            dr_iteracion_min = dr_iteracion;
        end        
        if( dt_iteracion < dt_iteracion_min)
            dt_iteracion_min = dt_iteracion;
        end        
    end
    
    iter_old = iter;
    t_iter(iter) = t; 
    iter = iter + 1;
end


% Mostrar valores finales:
numero_de_iteraciones = iter_old

pos_1_final = pos_1(iter_old,:)
pos_2_final = pos_2(iter_old,:)
vel_1_final = vel_1(iter_old,:)
vel_2_final = vel_2(iter_old,:)

Dif_v1 = norm(vel_1_final) - norm(v1_ini);
%Dif_v2 = norm(vel_2_final) - norm(v2_ini)

Sigma_1_final = 1 / sqrt(1 - norm(vel_1_final)^2);
Sigma_2_final = 1 / sqrt(1 - norm(vel_2_final)^2);

Sigma_1_inicial =  1 /sqrt(1 - norm(v1_ini)^2);
Sigma_2_inicial =  1 /sqrt(1 - norm(v2_ini)^2);

% Ignorando las masas, esta fórmula es la buena
Dif_Ener_cinet_1 = m1 * (Sigma_1_final - Sigma_1_inicial) 
Dif_Ener_cinet_2 = m2 * (Sigma_2_final  - Sigma_2_inicial)

% Mostramos las velocidades.
v_max_1 = v_max_1;
v_max_2 = v_max_2;
v_min_1 = v_min_1;
v_min_2 = v_min_2;


% Sólo válido para casos simétricos contra otro electrón
Sigma_v_min_1 = 1 / sqrt(1 - v_min_1^2);
E_min_1 = (Sigma_v_min_1 - 1) * m1;

Sigma_v_ini_1 = 1 / sqrt(1 - norm(v1_ini)^2);
E_ini_1 = (Sigma_v_ini_1 - 1) * m1;

Dif_Ec_max_1 = E_ini_1 - E_min_1;

[V_fin_1, A_fin_1] = Potencial_LW( (pos_1(iter_old,:)-pos_ret1), vel_ret1);
V_fin_1 = V_fin_1 * prod_q
A_fin_1 = A_fin_1 * prod_q;
Lagr_fin_1 = V_fin_1 - dot(A_fin_1, vel_1_final);

Dif_pot = V_fin_1 - V_ini_1
%Dif_lagr = Lagr_fin_1 - Lagr_ini_1


if consid_rad == 1
    %Var_Ener_Rad_1 = Var_Ener_Rad_1
    %Var_Ener_Rad_2 = Var_Ener_Rad_2
    Ener_Larmor_1 = Ener_Larmor_1 
    Ener_Larmor_2 = Ener_Larmor_2
end

L_1_ini = cross(pos_1_ini, v1_ini) 
L_1_final = cross(pos_1_final, vel_1_final)

L_2_ini = cross(pos_2_ini, v2_ini) 
L_2_final = cross(pos_2_final, vel_2_final)

dr_iteracion_min = dr_iteracion_min
dt_iteracion_min = dt_iteracion_min
