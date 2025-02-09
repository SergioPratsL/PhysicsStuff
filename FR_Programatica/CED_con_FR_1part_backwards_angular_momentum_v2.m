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

pos_ini = [1000000, 1, 0];

% Este es el parámetro más importante.
%v_ini = [0.5, 0, 0];      % Bien
%v_ini = [0.6, 0, 0];       % Bien
%v_ini = [0.7, 0, 0];      % Arreglado
%v_ini = [0.75, 0, 0];     
v_ini = [0.8, 0, 0];       % 0.9993
p_ang_ini = cross(pos_ini, v_ini);
 

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

factor_inversor = - 1;

dr_iteracion = d_ini / num_it_octava / 2;
% El tiempo va para atrás en este algoritmo!
dt_iteracion = factor_inversor * dr_iteracion / speed;

%Ener_Larmor1 = 0;

Ener_Larmor = 0;
p_ang_Larmor = [0,0,0];

pos = pos_ini;
pos_old = pos;
vel = v_ini;
vel_old = vel;
vel_old2 = vel;
mom = 1 / sqrt(1 - norm(vel)^2) * vel;
mom_old = mom;

iter = 2;       % Se incrementa al principio
t = 0;
% Simplemente es para forzar que se calcule la posición retardada a partir

speed = norm(vel(1,:));

Sigma = fGamma(vel(1,:));        % Así se calcula solo!
Sigma_old = Sigma;
mom_u = [Sigma, mom];
mom_u_old = mom_u;

Integral_DF_u_prop_old = [0,0,0,0];
acel_ext_prop_old = [0,0,0];
a_prop = [0,0,0];
a_prop_sqr = 0;
dv_prop_old = [0, 0, 0];
v_0 = [0,0,0];
revisar_dt_iteracion = 0;
dv_mom_old_cuadr = 0;

[factor_atenuacion, factor_inicial, ds_iteracion] = CalculaFactorAtenuacion(dt_iteracion, Sigma); 

% debug
fr_max = [0,0,0];
dist_min = 10000;

while dist < (d_ini_original + 1)
    % dt_iteracion es negativa por lo que el tiempo irá para atrás
    
    pos = pos_old + vel_old * dt_iteracion;  
    
    [E, B] = CampoInducido_sin_unidades(pos, v_0);    
    
    F_ext = FLorentz(vel_old, E, B);
    
    % La promiscuidad entre fuerza y aceleración se debe a que m=1
    a_u_ext = Sigma * [dot(F_ext, vel_old), F_ext];      
    a_u_ext_prop = Boost(a_u_ext, vel_old);
    acel_ext_prop = a_u_ext_prop(2:4);    
    
    if iter == 24000
        iter = iter;
    end
    
    % Espero que haciendo la rotación de las narices se resuelva este
    % entuerto! La afectación es mínima!
    if ( pos_ini(2) ~= 0 && dv_mom_old_cuadr > 0.000000001 )
        % Este paso tiene un gran coste computacional...
        rotMatrix = GetThomasRotMatrix( vel_old, vel_old2 );
        
        Integral_DF_u_prop_old(2:4) = Integral_DF_u_prop_old(2:4) * rotMatrix;
        % Esto influye mucho más que rotar la integral!
        acel_ext_prop_old = acel_ext_prop_old * rotMatrix;
    end
    
    Integral_DF_u_prop = Boost(Integral_DF_u_prop_old, dv_prop_old);
    
    dif_acel_ext_prop = (acel_ext_prop - acel_ext_prop_old) / ds_iteracion;
    
    % El (2/3) es porque la integral es (3/2)dF_ext + 1*a^2*a...
    integral_dF_prop = Integral_DF_u_prop(2:4) * factor_atenuacion + (dif_acel_ext_prop - (2/3) * a_prop_sqr * a_prop) * factor_inicial;

    a_prop = acel_ext_prop + (2/3) * integral_dF_prop;
    a_prop_sqr = a_prop(1)^2 + a_prop(2)^2 + a_prop(3)^2;
    
    % Este es el autovector de la aceleración propia, sabemos que en el SRI
    % origen la potencia es 0.
    a_prop_u = [0, a_prop];
    
    % Esto se pasará a la siguiente iteración para que haga boost, en el
    % SRI propio el término v·da = 0.
    Integral_DF_u_prop_old = [a_prop_sqr, integral_dF_prop];
    
    dv_mom_old = ds_iteracion * a_prop;
    dv_mom_old_cuadr = dv_mom_old(1)^2 + dv_mom_old(2)^2 + dv_mom_old(3)^2 ;
    if(dv_mom_old_cuadr > 0.000001)
        dv_Sigma = sqrt(1 + dv_mom_old_cuadr);
        dv_prop_old = dv_mom_old / dv_Sigma;
        
        % Tenemos que revisar el número de iteraciones porque no puede ser
        % que cambie tanto la velocidad entre iteraciones seguidas!
        revisar_dt_iteracion = 1;
    else
        dv_prop_old = dv_mom_old;
    end
    
    % Por fin se vuelve al SRI lab, se usa vel_old porque en este punto aún
    % no tenemos la nueva velocidad... (poner vel sólo sería engañarme
    % porque en este punto, vel = vel_old).
    a_u = Boost(a_prop_u, -vel_old);  % Probar transform velocidades

    % Se usa ds_iteracion porque aunque el intervalo es dt_iteracion, a_u
    % tiene un factor Sigma respecto a la fuerza y así se asnulan
    dif_p = a_u(2:4) * ds_iteracion;
    
    mom = mom_old + dif_p;
        
% v^2(1+p^2) = p^2  (m=1)
    mom_cuadr = mom(1)^2 + mom(2)^2 + mom(3)^2;
% E^2 = 1 + p^2 = Sigma^2   (m=1)
    Sigma = sqrt(1 + mom_cuadr);
    
    mom_u = [Sigma, mom];
    
    vel = mom / Sigma;
    
    % Aceleración del laboratorio, necesario para calcular Larmor.
    acel_lab = (vel - vel_old) / dt_iteracion;
    
    % Redundante, para comparar (coincide)
    %pot_rad = Larmor_covariante(vel, acel_lab, Sigma);
    %Ener_Larmor1 = Ener_Larmor1 + pot_rad * abs(dt_iteracion); 
    
    pot_rad_prop = Larmor_no_covariante(a_prop);
    % Siendo ortodoxos habría que usar Sigma * ds_iteracion, pero equivale
    % a dt_iteracion
    Ener_Larmor = Ener_Larmor + pot_rad_prop * abs(dt_iteracion);
    % fuerza_larmor = Sigma * pot_rad_prop * vel
    torque_Larmor = cross(pos, vel) * pot_rad_prop;
    p_ang_Larmor = p_ang_Larmor + torque_Larmor * abs(dt_iteracion);
    
    
    % Para ahorrar cálculos no evaluaré el factor de iteracción cada
    % vez sino sólo cuando la Sigma hayacambiado más de una milésima
    if abs(Sigma - Sigma_old) > 0.0001
        [factor_atenuacion, factor_inicial, ds_iteracion] = CalculaFactorAtenuacion(dt_iteracion, Sigma); 
        Sigma_old = Sigma;
    end
    
% Por último evaluar si hay que cambiar de octava.
    dist = norm(pos);
    
    if dist < dist_min 
        dist_min = dist;
    end
    
    if dist <= (d_ini / 2) || dist >= (d_ini * 2) || revisar_dt_iteracion == 1
       
        if revisar_dt_iteracion == 1
            % Reducir las iteraciones a la mitad las veces que haga falta
            dr_iteracion = dr_iteracion / 2;
            revisar_dt_iteracion = 0;
        else
            dr_iteracion_test = dist / num_it_octava;
            % Prudencia no quiero forzar un error, pero esto se ha de
            % quitar cuando las aceleraciones vuelven a ser débiles.
            if dr_iteracion_test > (2 * dr_iteracion) && norm(a_u(2:4)) > 0.01   
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
    vel_old2 = vel_old;
    vel_old = vel;
    mom_old = mom;
    acel_ext_prop_old = acel_ext_prop;
    mom_u_old = mom_u;
    
    iter = iter + 1;
end

% Mostrar valores finales:
numero_de_iteraciones = iter-1
iter;
pos_final = pos_old
vel_final = vel_old

p_ang_fin = cross(pos_final, vel_final)

dist_min = dist_min

Dif_v1 = norm(vel_final) - norm(v_ini);

Sigma_final = 1 / sqrt(1 - norm(vel_final)^2);
Sigma_inicial =  1 /sqrt(1 - norm(v_ini)^2);
Dif_Ener_cinet_1 = (Sigma_final - Sigma_inicial) 

% Error para v=0.4 (sin FR)
%Dif_Ener_cinet_1_real = Dif_Ener_cinet_1 - 1.5322e-05;

%Ener_Larmor1 = Ener_Larmor1 
Ener_Larmor = Ener_Larmor
p_ang_Larmor = p_ang_Larmor

ratio_conservacion_mom_ang = norm(p_ang_fin)  / norm(p_ang_ini + p_ang_Larmor)

L_1_ini = cross(pos_ini, v_ini);
L_1_final = cross(pos_final, vel_final);

dr_iteracion_min = dr_iteracion_min;
dt_iteracion_min = dt_iteracion_min;

Ratio_Dif_Ener_Larmor = Dif_Ener_cinet_1 / Ener_Larmor
%Ratio_Dif_Ener_Larmor_real = Dif_Ener_cinet_1_real / Ener_Larmor

%fr_max = fr_max
