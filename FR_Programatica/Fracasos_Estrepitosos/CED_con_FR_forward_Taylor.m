% Script para evaluar la fuerza de reacción sobre una carga que va a gran
% velocidad acercándose a otra que está estática y no acelera (su masa es
% infinita). La carga 2 está quieta en el centro.

% Este script se basa en el que iba hacia atras en el tiempo 
% pero esta vez usaré Taylor para ver si con eso podemos obtener buena
% aproximación de la derivada de la aceleraciónn y no tener que ir hacia
% atrás en el tiempo.

% El resultado final es que la aproximación de Taylor no ayuda nada o tal
% vez haga falta una cantidad muy grande de términos, pero yo no he podido
% lograr que esto vaya bien.

% No hay muchas constantes en este matlab, la permitividad y permeabilidad
% están fuera, la carga y la masa también y la velocidad de la luz tb.
% en es

% m1 = 1;
% q1 = 1;
% c = 1
% 4pi_epsilon = 1

clear 

pos_ini = [-1000000, 0, 0];

% Este es el parámetro más importante.
%v_ini = [0.4, 0, 0];      
%v_ini = [0.5, 0, 0];      
%v_ini = [0.6, 0, 0];      
v_ini = [0.7, 0, 0];      
%v_ini = [0.75, 0, 0];    
%v_ini = [0.8, 0, 0];     

d_ini_original = norm(pos_ini);   % Es la distancia a la que empezó la prueba...
d_ini = d_ini_original;             % Controla el cambio de octava
dist = d_ini;

% Para saber cuánto se mueven por unidad de tiempo, en un primer momento
speed = norm(v_ini);

% Iteraciones por octava (es decir, hasta que se divida o multiplique
% por dos la distancia inicial
%num_it_octava = 2000;
num_it_octava = 4000;

dr_iteracion = d_ini / num_it_octava / 2;
dt_iteracion = dr_iteracion / speed;

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

Sigma = fGamma(vel(1,:));      
ds_iteracion = dt_iteracion / Sigma;
ds_iteracion_old = ds_iteracion;

a_ext_prop_old = [0,0,0];
a_prop = [0,0,0];
a_prop_old = [0,0,0];
v_0 = [0,0,0];
revisar_dt_iteracion = 0;

da_prop = [0,0,0];
d2a_prop = [0,0,0];
d3a_prop = [0,0,0];
d4a_prop = [0,0,0];
d5a_prop = [0,0,0];

da_prop_old = [0,0,0];
d2a_prop_old = [0,0,0];
d3a_prop_old = [0,0,0];
d4a_prop_old = [0,0,0];
d5a_prop_old = [0,0,0];

% debug
dist_min = 10000;

% Debug
dr_iteracion_min = 10000000;
dt_iteracion_min = 10000000;

while dist < (d_ini_original + 1)
    
    pos = pos_old + vel_old * dt_iteracion;  
    
    [E, B] = CampoInducido_sin_unidades(pos, v_0);    
    
    F_ext = FLorentz(vel_old, E, B);
    
    if iter == 16170 || iter == 16185
        iter = iter;
    end        
        
    % La promiscuidad entre fuerza y aceleración se debe a que m=1
    a_u_ext = Sigma * [dot(F_ext, vel_old), F_ext];      
    a_u_ext_prop = Boost(a_u_ext, vel_old);
    a_ext_prop = a_u_ext_prop(2:4);    
    
%     % Taylor y no los resultados darán la derivada de la aceleración dudo mucho que ayude.  
%     if(iter == 2)
%         da_prop = (a_prop - a_prop_old) / ds_iteracion_old;
%     else
%         da_prop = AproximacionTaylor5( da_prop, d2a_prop, d3a_prop, d4a_prop, d5a_prop, ds_iteracion_old);
%     end

    % último intento (3º) calcular a_prop_est con Tylor
    a_prop_est = AproximacionTaylor5( a_prop, da_prop, d2a_prop, d3a_prop, d4a_prop, ds_iteracion);
    da_prop = (a_prop_est - a_prop_old) / ds_iteracion;
    

    a_prop = a_ext_prop + (2/3) * da_prop;
    
    % Hacemos una aproximación de Taylor de la derivada de la aceleración
    %a_prop = a_ext_prop + (2/3) * AproximacionTaylor5( da_prop, d2a_prop, d3a_prop, d4a_prop, d5a_prop, ds_iteracion);
    %da_prop = (a_prop - a_prop_old) / ds_iteracion;
    
    if(iter > 3)
        d2a_prop = (da_prop - da_prop_old) / ds_iteracion_old;
    end
    
    if(iter > 4)
        d3a_prop = (d2a_prop - d2a_prop_old) / ds_iteracion_old;
    end    
    
    if(iter > 5)
        d4a_prop = (d3a_prop - d3a_prop_old) / ds_iteracion_old;
    end            
    
    if(iter > 6)
        d5a_prop = (d4a_prop - d4a_prop_old) / ds_iteracion_old;
    end                
    
    ds_iteracion_old = ds_iteracion;
    
    % Este es el autovector de la aceleración propia, sabemos que en el SRI
    % origen la potencia es 0.
    a_prop_u = [0, a_prop];
    
    dv_mom_old = ds_iteracion * a_prop;
    dv_mom_old_cuadr = dv_mom_old(1)^2 + dv_mom_old(2)^2 + dv_mom_old(3)^2 ;
    % Tenemos que revisar el número de iteraciones porque no puede ser
    % que cambie tanto la velocidad entre iteraciones seguidas!
    if(dv_mom_old_cuadr > 0.000001)        
        revisar_dt_iteracion = 1;
    end
    
    % Por fin se vuelve al SRI lab, se usa vel_old porque en este punto aún
    % no tenemos la nueva velocidad... (poner vel sólo sería engañarme
    % porque en este punto, vel = vel_old).
    a_u = Boost(a_prop_u, -vel_old);
    
    % Se usa ds_iteracion porque aunque el intervalo es dt_iteracion, a_u
    % tiene un factor Sigma respecto a la fuerza y así se asnulan
    dif_p = a_u(2:4) * ds_iteracion;
    
    mom = mom + dif_p;
        
% v^2(1+p^2) = p^2  (m=1)
    mom_cuadr = mom(1)^2 + mom(2)^2 + mom(3)^2;
    
% E^2 = 1 + p^2 = Sigma^2   (m=1)
    Sigma = sqrt(1 + mom_cuadr);
    
    vel = mom / Sigma;
    
    % Aceleración del laboratorio, necesario para calcular Larmor.
    acel_lab = (vel - vel_old) / dt_iteracion;
    
    pot_rad = Larmor_covariante(vel, acel_lab, Sigma);
    Ener_Larmor = Ener_Larmor + pot_rad * dt_iteracion; 
        
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
        dt_iteracion = dr_iteracion / speed;
        ds_iteracion = dt_iteracion / Sigma;        
        
        if( dr_iteracion < dr_iteracion_min)
            dr_iteracion_min = dr_iteracion;
        end        
        if( abs(dt_iteracion) < abs(dt_iteracion_min))
            dt_iteracion_min = abs(dt_iteracion);
        end
    end
    
    pos_old = pos;
    vel_old = vel;
    
    a_prop_old = a_prop;
    da_prop_old = da_prop;
    d2a_prop_old = d2a_prop;
    d3a_prop_old = d3a_prop;
    d4a_prop_old = d4a_prop;
    d5a_prop_old = d5a_prop;
    
    if iter >= 21000 && (mod(iter, 1000) == 0)
        iter = iter;
    end    
    
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

Ratio_Dif_Ener_Larmor = - Dif_Ener_cinet_1 / Ener_Larmor

