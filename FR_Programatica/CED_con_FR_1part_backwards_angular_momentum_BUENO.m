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
%v_ini = [0.5, 0, 0];      % bien  
%v_ini = [0.6, 0, 0];       % Lo peta!
%v_ini = [0.7, 0, 0];     % bien
%v_ini = [0.75, 0, 0];     
v_ini = [0.8, 0, 0];      % Lo peta

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
vel = v_ini;
% vel_old es usada para medir la aceleración en el laboratorio para Larmor
vel_old  = vel;
mom = 1 / sqrt(1 - norm(vel)^2) * vel;

iter = 2;       % Se incrementa al principio
t = 0;
% Simplemente es para forzar que se calcule la posición retardada a partir

speed = norm(vel(1,:));

Sigma = fGamma(vel(1,:));        % Así se calcula solo!
Sigma_old = Sigma;
Sigma_vel = Sigma * vel;

Integral_lab = [0,0,0];
v_0 = [0,0,0];
revisar_dt_iteracion = 0;

% Fuerza total de las dos últimas iteraciones para sacar ds_a que
% lamentablemente necesitamos en este caso al no fiarnos de los múltiples
% boosts.
F_u_old = [0,0,0,0];
F_u_old_old = [0,0,0,0];

F_ext_old = [0,0,0];

[factor_atenuacion, factor_inicial, ds_iteracion] = CalculaFactorAtenuacion(dt_iteracion, Sigma); 

dist_min = 10000;

while dist < (d_ini_original + 1)
    % dt_iteracion es negativa por lo que el tiempo irá para atrás
    
    % Es la velocidad obtenida de la anterior iteracion
    pos = pos + vel * dt_iteracion;  
    
    [E, B] = CampoInducido_sin_unidades(pos, v_0);    
    
    F_ext = FLorentz(vel, E, B);
    
    ds_F_ext = (F_ext - F_ext_old) / ds_iteracion;
    
    termino_F_ext = Sigma*dot(vel, F_u_old(2:4))*F_ext + Sigma*ds_F_ext;
   
    a_cuadr = Sigma^2 * ( F_u_old(1)^2 - (F_u_old(2)^2 + F_u_old(3)^2 + F_u_old(4)^2));
    
    % Este cabron me va a hundir
    % ds_a_u = Sigma * (F_u_old - F_u_old_old) / ds_iteracion;
    ds_F_u = (F_u_old - F_u_old_old) / ds_iteracion;
    ds_a_u = Sigma*dot(vel, F_u_old(2:4))*F_u_old + Sigma * ds_F_u;
    
    % Este término es la fuente de todos los problemas   
    prod_a_u_ds_a_u = Sigma * ( F_u_old(1)*ds_a_u(1) - dot(F_u_old(2:4), ds_a_u(2:4)) );
    
    termino_hijo_de_la_gran_puta = Sigma * a_cuadr * F_u_old(2:4) + 2 * prod_a_u_ds_a_u  * Sigma_vel;
    
    % El (2/3) es porque la integral es (3/2)dF_ext + 1*a^2*a...
    Integral_lab = Integral_lab * factor_atenuacion + (termino_F_ext + (2/3) * termino_hijo_de_la_gran_puta) * factor_inicial;

    F = F_ext + (2/3) * (Integral_lab + a_cuadr * Sigma_vel) / Sigma;

    % Se usa ds_iteracion porque aunque el intervalo es dt_iteracion, a_u
    % tiene un factor Sigma respecto a la fuerza y así se asnulan
    dif_p = F * dt_iteracion;
    
    
    dif_p_cuadr = dif_p(1)^2 + dif_p(2)^2 + dif_p(3)^2 ;
    if(dif_p_cuadr > 0.000001)
        revisar_dt_iteracion = 1;
    end

    % Y por fin actualizamos el momento y la velocidad :)
    mom = mom + dif_p;
    
% v^2(1+p^2) = p^2  (m=1)
    mom_cuadr = mom(1)^2 + mom(2)^2 + mom(3)^2;
% E^2 = 1 + p^2 = Sigma^2   (m=1)
    Sigma = sqrt(1 + mom_cuadr);
    
    vel = mom / Sigma;
    
    % Ya que podemos, calculamos la F_u con la velocidad actualizada.
    F_u = [dot(F, vel), F];
    
    % Aceleración del laboratorio, necesario para calcular Larmor.
    acel_lab = (vel - vel_old) / dt_iteracion;
    
    pot_rad = Larmor_covariante(vel, acel_lab, Sigma);
    Ener_Larmor = Ener_Larmor + pot_rad * abs(dt_iteracion); 

    % La potencia radiada en el Lab es la misma que en el SRI propio ya que
    % el sigma del Boost de [E,0,0,0] se compensa con la dilatacion
    % del tiempo, es boost de la partícula a lab, por tanto -vel
    torque_Larmor = cross(pos, -vel) * pot_rad;
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
            if dr_iteracion_test > (2 * dr_iteracion) && norm(F) > 0.01   
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
    
    vel_old = vel;
    F_ext_old = F_ext;
    F_u_old_old = F_u_old;
    F_u_old = F_u;
    Sigma_vel = Sigma * vel;
    
    iter = iter + 1;
end

% Mostrar valores finales:
numero_de_iteraciones = iter-1
iter;
pos_final = pos
vel_final = vel

p_ang_fin = cross(pos_final, vel_final)

dist_min = dist_min

Dif_v1 = norm(vel_final) - norm(v_ini);

Sigma_final = 1 / sqrt(1 - norm(vel_final)^2);
Sigma_inicial =  1 /sqrt(1 - norm(v_ini)^2);


% Teniendo en cuenta que el el proceso va hacia atrás en el tiempo.
% Espero que siempre sea una cantidad positiva!!
Ec_Perdida = (Sigma_final - Sigma_inicial) 

%Ener_Larmor1 = Ener_Larmor1 
Ener_Larmor = Ener_Larmor
p_ang_Larmor = p_ang_Larmor

ratio_conservacion_mom_ang = norm(p_ang_Larmor)  / norm(p_ang_ini - p_ang_fin)

L_1_ini = cross(pos_ini, v_ini)
L_1_final = cross(pos_final, vel_final)

dr_iteracion_min = dr_iteracion_min;
dt_iteracion_min = dt_iteracion_min;

Ratio_EcPerdida_ELarmor =  Ec_Perdida / Ener_Larmor

