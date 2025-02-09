

iteraciones = 40000;

% Inicializo las variables
tiempo = zeros(1, iteraciones);

% Posiciones
r1 = zeros(3, iteraciones);
r2 = zeros(3, iteraciones);

% Momento inercial
p1 = zeros(3, iteraciones);
p2 = zeros(3, iteraciones);

% Velocidades
v1 = zeros(3, iteraciones);
v2 = zeros(3, iteraciones);

%Aceleraciones
a1 = zeros(3, iteraciones);
a2 = zeros(3, iteraciones); 

% 4-Potenciales, obviamente solo hay uno pues es compartido.
Potencial = zeros(4, iteraciones);


R = [0 0 0];


% En este caso las particulas van a chocar frontalmente (no deberian llegar
% a chocar pero bueno).
% Elegimos como SRI de referencia aquel en el cual las particulas se mueven
% a igual velocidad pero en sentido contrario y medimos desde la posicion
% en la que estan equidistantes.
% Tiene la ventaja de ser un escenario simetrico en el que lo que se evalua
% para una carga servia tambien para la otra

% Todo sera en unidades SRI para joder :P

% Velocidad de la luz m/s
c = 2.9979 * 10^8;

% Carga del electron (Coulombs
e_q = 1.6022 * 10^-19;

% Masa del electron (Kg)
e_m = 9.1094 * 10^-31;

% Constante de permitividad (Farads por metro)
E_0 = 8.8542 * 10^-12;

% Constante de permeabilidad (Newtons / Ampere^2)
M_0 = 4 * pi() * 10^-7


% Constante de aceleracion propia
% Esta constante daria la aceleracion en m/s a una distancia de 1m.
C_F = 1/(4*pi*E_0) * e_q^2 / e_m;

% Y esta es la aceleracion en velocidad de la luz por segundo
C_F_c = C_F / c;



% Las particula se mueven a velocidades relativistas, si no tiene gracia
v_ini = 0.2 * c;  % Siempre en positivo

% La distancia inicial sera un metro
r_ini = 0.5;      % Siempre en positivo

% Fraccion de tiempo que avanzara cada vez
% En cada momento esto sera proporcional a la distancia entre los
% electrones
fraccion_base = 1 / 100;


% Energia cinetica de cada de las particulas
Sigma = fGamma(v_ini/c);
Ec = e_m * c^2 * (Sigma - 1);

% La energia cinetica total sera la suma de las dos, por tanto el doble
Ec_tot = 2 * Ec;

% Para hacerme una idea de a que distancia se cancelaria el potencial lo
% uso Coulomb a saco

R_min_aprox = e_q^2 / (4*pi* E_0 * Ec_tot);

% Para hacerme una idea: haria falta llegar a 6.8331e-014 m para cancelar
% toda la energia, esto es menos que la longitud de compton del electron:
%  2.4263102389(16)×10?12

% En todo caso en este escenario estamos exentos de esas disquisiciones.


salir = 0;

% Nos moveremos en la X por esta vez
iter = 1;
t = 0;
tiempo(1, iter) = t;

r1(1, iter) = -r_ini;
r2(1, iter) = r_ini;

v1(1, iter) = v_ini;
v2(1, iter) = -v_ini;

p1(1, iter) = v_ini * fGamma(v_ini/c) * e_m;
p2(1, iter) = -v_ini * fGamma(v_ini/c) * e_m;

% Energia y momento potenciales inicial, usamos el de L-W a falta de algo mejor:
% Potencial(1,iter) = 1/(4*pi*E_0) * e_q^2 /(2*r_ini) / (1-(v_ini/c));
% Potencial(2,iter) = 1/(4*pi*E_0) * e_q^2 /(2*r_ini) / (1-(v_ini/c)) * (v_ini/c);
% Como solo interesan las diferencias derivadas del proceso, los dejare a cero.



% Las formulas de este loop solo son validas para caso de choque frontal

% Inicializo esta variable cuya gestion va a darme muchos problemas pues
% representa el momento desde el que se recibe el campo que evaluas en la
% otra particula
iter_ret = 1;

while salir == 0
   
% Vector que va de "2" a "1"        
    r1_act = r1(:, iter);
    r2_act = r2(:, iter);
    
    R_act = r1_act - r2_act;    
    dist_act = norm(R_act);    

    
% PUNTO CLAVE: Hay que sacar la distancia pasada (esto seria bueno
% encapsularlo en una rutina) 

% Si el instante pasado no lo tenemos, hacemos una aproximacion a velocidad
% constante
    if t - dist_act / (c - v_ini ) < 0
        iter_ret = 1;
        dist_ret = dist_act / (1 - v_ini / c);
        v_ret = v_ini;
        a_ret = 0;      % Poco relevante.
        
    else
% Hay que encontrar los intervalos en los que se cumple que el tiempo
% retardado

% Primero ver en que sentido debo moverme para localizar la iteracion buena
% (generalmente sera hacia adelantan
        
% La variable que sale del LOOP es iter_ret
        buscando = 'X'; 
        primero = 'X';
        while buscando == 'X'
        
            lapso1 = tiempo(iter) - tiempo(iter_ret);
        
% El vector retardado contiene la distancia desde el emisor en el tiempo
% retardado hasta el receptor en el tiempo actual        
            r2_ret1 =r2(:, iter_ret);        
            R_ret1 = r1_act - r2_ret1;
            dist_ret1 = norm( R_ret);
            t_vuelo1 = dist_ret / c;
            
            dif_lapso_t_vuelo_1 = lapso1 - t_vuelo1;
        
% La diferencia de tiempo entre los dos momentos cogidos
% es menor que la distancia (con c normalizado), hay que ir hacia atras para
% lograr que la diferencia de instantes iguale a la distancia
            if dif_lapso_t_vuelo_1 <= 0
        
                if iter_ret > 1
                    iter_ret = iter_ret - 1;
                else
                    buscando = '';
                    break
                end

% La diferencia de tiempo entre los dos momentos cogidos
% es mayor que la distancia (con c normalizado), hay que ir hacia atras para
% lograr que la diferencia de instantes iguale a la distancia
            elseif dif_lapso_t_vuelo_1 > 0

                if iter_ret < iter
                    iter_ret = iter_ret + 1;
                else
                    buscando = '';
                    break
                end                
                  
            end    

% Vemos si ya hemos encontrado el segmento bueno en el cual salio la luz
% que ahora llega
            if primero == 'X'
                primero = '';
            else
                if ( dif_lapso_t_vuelo_old <= 0 && dif_lapso_t_vuelo_1 >= 0 ) || ( dif_lapso_t_vuelo_old >= 0 && dif_lapso_t_vuelo_1 <= 0 )

% Nos quedamos con el valor que cuya distancia recorrida se aproxime mas a
% la diferencia de tiempo. Con iter_ret o con el valor de iter_ret que
% estabamos mirando antes
                    if abs(dif_lapso_t_vuelo_old) < abs(dif_lapso_t_vuelo_1)
                        iter_ret = iter_old; % La iteracion que se ajusta mas es la anterior
                    end
                    
                    buscando = '';
                    break   % salimos del while
                end
            end
            
            dif_lapso_t_vuelo_old = dif_lapso_t_vuelo_1;
% El valor anterior que estabamos mirando de este campo
            iter_old = iter_ret;        
            
            if iter_ret == 1 || iter_ret == iter
                buscando = '';                
                break       % Sale del LOOP
            end            
            
        end     % while
    end 
  
% Ahora si, obtenemos el tiempo y posicion retardadas    
    r2_ret = r2(:, iter_ret);        
    R_ret = r1_act - r2_ret;
    dist_ret = norm( R_ret);   
    v2_ret = v2(:,iter_ret);
    a2_ret = a2(:,iter_ret); % realmente todavia no nos hace falta

%~PUNTO CLAVE: Hay que sacar la distancia pasad    
    
% Hacemos un seguimiento todavia mas exhaustivo en los momentos en que
% ambas particulas estan muy cerca
    if dist_act > 10^-7
        CoefN = 1;
    elseif dist_act > 10^-9
        CoefN = 2;
    elseif dist_act > 10^-11
        CoefN = 3;
    else
        CoefN = 4;
    end
    
    dif_t = ( fraccion_base / CoefN )  * (dist_ret / c);
    t = t + dif_t;
    
% Por el momento coger la distancia anterior sin considerar lo que se pueda
% acelerar en el intervalo
    v = v1(:, iter);
    dif_pos = v * dif_t;
    
% Campo electrico sobre la otra particula (podemos descartar el magnetico)
    [E, B] = CampoInducido(R_ret, v2_ret, e_q);
    
% Aceleracion en metros por segundo, como no aplica el magnetismo en el escenario de acercamiento frontal no lo pondre    
% hay que vigilar, puesto que al ser otro SRI no podemos hablar de
% aceleracion propiamente dicha
    dp1_dt = E * e_q;

% Variacion de momento    
    dp = dp1_dt * dif_t;
    
    p1_act = p1(:,iter) + dp;
    p1_massless_norm = p1_act / e_m / c;

    v1_modul = sqrt( norm(p1_massless_norm)^2 / (1 + norm(p1_massless_norm)^2) ) * c;
    
% Una vez tenemos el modulo, la direccion se recupera de p1_act    
    p1_unitario = p1_act / norm(p1_act);
    
    v1_act = v1_modul * p1_unitario;
    
% Por fin la aceleracion, imprescindible para saber cuanta energia se radia en el infinito
% Es la aceleracion vista desde el laboratorio, no la aceleracion propia
    a1_act = ( v1_act - v1(iter) ) / dif_t;
    

%  El potencial vector sera cero por simetria    
    dE = 2 * dot(dp, v1(:,iter));
    V_act = Potencial(1,iter) + dE;
       
    
% Escribimos los valores    
    iter = iter + 1;
    
    tiempo(1,iter) = t;
    
    r1(:,iter) = r1_act + dif_pos;
    r2(:,iter) = r2_act - dif_pos;
    
    p1(:, iter) = p1_act + dp;
    p2(:, iter) = -p1_act - dp;  % Se deduce.
    
    v1(:,iter) = v1_act;
    v2(:,iter) = -v1_act;       % Se deduce
    
    a1(:,iter) = a1_act;
    a2(:,iter) = -a1_act;
        
    Potencial(1,iter) = V_act;    

%    Condicion de salida (que las particulas hayan rebotado y superen la
%    distancia inicial o bien que nos quedemos sin array (algo ha ido mal :S)
    if iter >= iteraciones || norm( r1_act ) > r_ini
        salir = 1;
    end
    
end

% Mostramos en que iteracion se quedo
iter


% Mostramos el potencial que queda al final
Potencial(0,iter)

% Es algo obvio, pero esta bien ver que dicho potencial valga lo mismo que
% la variacion de energia cinetica
Sigma_f  = fGamma( v1(:,iter)/c );

Ec_fin = 2 * e_m * c^2 * (Sigma_f - 1)

Dif_Ec = Ec_fin - Ec_tot






