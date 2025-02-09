% Script para evaluar la fuerza de reacción sobre una carga que va a gran
% velocidad acercándose a otra que está estática y no acelera (su masa es
% infinita). La carga 2 está quieta en el centro.

% Este script va hacia atras en el tiempo para así depender de valores de
% la derivada de la aceleración que ya hemos recorrido

% No hay muchas constantes en este matlab, la permitividad y permeabilidad
% están fuera, la carga y la masa también y la velocidad de la luz tb.
% en es

% m1 = 1;
% q1 = 1;
% c = 1
% 4pi_epsilon = 1

clear

pos_ini = [1000000, 0, 0];

% Este es el parámetro más importante.
v_ini = [0.4, 0, 0];      % ds_a =   1.0e+99 *    -8.4593         0         0
% nada más quie decir.

d_ini_original = norm(pos_ini);   % Es la distancia a la que empezó la prueba...
d_ini = d_ini_original;             % Controla el cambio de octava
dist = d_ini;

% Prueba de momento angular...
%pos_ini = [1000000, 1, 0];

% debug
dist_old = d_ini_original;

% Para saber cuánto se mueven por unidad de tiempo, en un primer momento
speed = norm(v_ini);

% Debug
dr_iteracion_min = 10000000;
dt_iteracion_min = 10000000;

% Iteraciones por octava (es decir, hasta que se divida o multiplique
% por dos la distancia inicial
num_it_octava = 2000;                
%num_it_octava = 4000;

% Estimar el tamaño de los vectores para no tener que recortarlos (memoria
% no debería faltar para un algoritmo como este)
num_octavas_estimado = round(log2(d_ini_original)) + 20;      

tam_vector = num_octavas_estimado * num_it_octava * 2;

factor_inversor = - 1;

dr_iteracion = d_ini / num_it_octava / 2;
% El tiempo va para atrás en este algoritmo!
dt_iteracion = factor_inversor * dr_iteracion / speed;

Ener_Larmor = 0;

pos = pos_ini;
pos_old = pos;
vel = v_ini;
vel_old = vel;
mom = 1 / sqrt(1 - norm(vel)^2) * vel;

iter = 2;       % Se incrementa al principio
t = 0;
% Simplemente es para forzar que se calcule la posición retardada a partir

speed = norm(vel(1,:));

Sigma = fGamma(vel(1,:));        % Así se calcula solo!
Sigma_old = Sigma;

acel = [0,0,0];
acel_old = acel;
acel_old_old = acel;
dv_old = [0, 0, 0];
v_0 = [0,0,0];
revisar_dt_iteracion = 0;

[factor_atenuacion, factor_inicial, ds_iteracion] = CalculaFactorAtenuacion(dt_iteracion, Sigma); 
ds_iteracion_old = ds_iteracion;

% debug
fr_max = [0,0,0];
dist_min = 10000;

while dist < (d_ini_original + 1)
    % dt_iteracion es negativa por lo que el tiempo irá para atrás
    
    pos = pos_old + vel * dt_iteracion;  
    
    [E, B] = CampoInducido_sin_unidades(pos, v_0);    
    
    F_ext = FLorentz(vel, E, B);
     
    if iter == 12070
        iter = iter;
    end
    
    % Hace falta Sigma para mostrar las componentes de a_u, y tb
    % ds_iteracion
    ds_a = Sigma * (acel_old - acel_old_old) / ds_iteracion_old;

    % 1/Sigma^2 es el efecto de restar el término de la potencial...
    % pero se compensa con que a_u tiene un factor Sigma...
    a_cuadr = (acel_old(1)^2 + acel_old(2)^2 + acel_old(3)^2);
    
    acel = F_ext + (2/3) * (ds_a - a_cuadr * vel);
    
    dif_p = acel * dt_iteracion;
    
    % Controlar si entre dos instantes la diference de fuerza es demasiado
    % grande
    dif_p_cuadr = dif_p(1)^2 + dif_p(2)^2 + dif_p(3)^2 ;
    if(dif_p_cuadr > 0.000001)
        revisar_dt_iteracion = 1;
    end
    
    % Prescindir de mom_old...
    mom = mom + dif_p;
        
% v^2(1+p^2) = p^2  (m=1)
    mom_cuadr = mom(1)^2 + mom(2)^2 + mom(3)^2;
    
% E^2 = 1 + p^2 = Sigma^2   (m=1)
    Sigma = sqrt(1 + mom_cuadr);
    
    mom_u = [Sigma, mom];
    mom_u_old = mom_u;
    
    % vel_old se actualiza antes de calcular la nueva velocidad.
    vel_old = vel;
    vel = mom / Sigma;

    % 4 vector de la aceleración necesasrio para la fuerza de reacción
    acel_u = Sigma * (mom_u - mom_u_old) / dt_iteracion;    
    
    % Aceleración del laboratorio, necesario para calcular Larmor.
    acel = (vel - vel_old) / dt_iteracion;
    
    pot_rad = Larmor_covariante(vel, acel, Sigma);
    %pot_rad = Larmor_no_covariante(acel);
    Ener_Larmor = Ener_Larmor + pot_rad * abs(dt_iteracion); 
    
    % Para ahorrar cálculos no evaluaré el factor de iteracción cada
    % vez sino sólo cuando la Sigma hay acambiado más de una milésima
    if abs(Sigma - Sigma_old) > 0.0003
        [factor_atenuacion, factor_inicial, ds_iteracion] = CalculaFactorAtenuacion(dt_iteracion, Sigma); 
        Sigma_old = Sigma;
    end
    
% Por último evaluar si hay que cambiar de octava.
    dist = norm(pos);
    
    if dist < dist_min 
        dist_min = dist;
    end
    
    % Para sacar dt_a con ls tiempo anteriores.
    ds_iteracion_old = ds_iteracion;
    
    if dist <= (d_ini / 2) || dist >= (d_ini * 2) || revisar_dt_iteracion == 1
       
        if revisar_dt_iteracion == 1
            % Reducir las iteraciones a la mitad las veces que haga falta
            dr_iteracion = dr_iteracion / 2;
            revisar_dt_iteracion = 0;
        else
            dr_iteracion_test = dist / num_it_octava;
            % Prudencia no quiero forzar un error, pero esto se ha de
            % quitar cuando las aceleraciones vuelven a ser débiles.
            if dr_iteracion_test > (2 * dr_iteracion) && norm(acel) > 0.01   
                dr_iteracion = 2 * dr_iteracion;
                d_ini = 2 * d_ini;
            else
                dr_iteracion = dr_iteracion_test;
                d_ini = dist;
            end
        end
        
        speed = norm(vel);
        dt_iteracion = factor_inversor * dr_iteracion / speed;
        
        [factor_atenuacion, factor_inicial, ds_iteracion] = CalculaFactorAtenuacion(dt_iteracion, Sigma); 
        
        if( dr_iteracion < dr_iteracion_min)
            dr_iteracion_min = dr_iteracion;
        end        
        if( abs(dt_iteracion) < abs(dt_iteracion_min))
            dt_iteracion_min = abs(dt_iteracion);
        end
    end
    
    pos_old = pos;
    acel_old_old = acel_old;
    acel_old = acel;   
    mom_u_old = mom_u;
    
    iter = iter + 1;
end

% Mostrar valores finales:
numero_de_iteraciones = iter-1
iter;
pos_final = pos_old
vel_final = vel_old

dist_min = dist_min

Dif_v1 = norm(vel_final) - norm(v_ini);

Sigma_final = 1 / sqrt(1 - norm(vel_final)^2);
Sigma_inicial =  1 /sqrt(1 - norm(v_ini)^2);
Dif_Ener_cinet_1 = (Sigma_final - Sigma_inicial) 

% Error para v=0.4 (sin FR)
%Dif_Ener_cinet_1_real = Dif_Ener_cinet_1 - 1.5322e-05;

Ener_Larmor = Ener_Larmor 

L_1_ini = cross(pos_ini, v_ini);
L_1_final = cross(pos_final, vel_final);

dr_iteracion_min = dr_iteracion_min;
dt_iteracion_min = dt_iteracion_min;

Ratio_Dif_Ener_Larmor = Dif_Ener_cinet_1 / Ener_Larmor
%Ratio_Dif_Ener_Larmor_real = Dif_Ener_cinet_1_real / Ener_Larmor

%fr_max = fr_max
