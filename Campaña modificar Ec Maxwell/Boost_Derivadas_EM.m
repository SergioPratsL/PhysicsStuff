function [E_uv_new, B_uv_new ] = Boost_Derivadas_EM(E_uv_ini, B_uv_ini, v_boost )
% Función para transformar las derivas del campo EM en un punto a otro SRI
% v_boost es un 3-vector con la velocidad normalizada
% De momento sólo está implementada para hacer boosts en dirección X!
% v_boost normalizado, c=1.

Sigma = 1 / sqrt(1-norm(v_boost)^2);

% Lo primero es obtener las derivadas de cada campo en cada dirección
E_dt_ini = E_uv_ini(1, 2:4); 
E_dx_ini = E_uv_ini(2, 2:4);
E_dy_ini = E_uv_ini(3, 2:4);
E_dz_ini = E_uv_ini(4, 2:4);

B_dt_ini = B_uv_ini(1, 2:4); 
B_dx_ini = B_uv_ini(2, 2:4);
B_dy_ini = B_uv_ini(3, 2:4);
B_dz_ini = B_uv_ini(4, 2:4);


% Siguiente sacar campos con en el nuevo SRI pero sin transformar las 
% coordenadas, me basaré en matrices más que en el viejo estilo, la línea
% compentada sólo sirve para recordar cual era ese viejo estilo

%[E_lab_dt_ori, B_lab_dt_ori] = Boost_EM(E_dt, B_dt, v_boost);

F_uv_dt_ini = TensorEM(E_dt_ini, B_dt_ini);
F_uv_dx_ini = TensorEM(E_dx_ini, B_dx_ini);
F_uv_dy_ini = TensorEM(E_dy_ini, B_dy_ini);
F_uv_dz_ini = TensorEM(E_dz_ini, B_dz_ini);

T_boost = Tensor_boosts(v_boost);   

% Transformo el vector al nuevo SRI
F_uv_new_dt_ini = (T_boost * F_uv_dt_ini) * T_boost' ; 
F_uv_new_dx_ini = (T_boost * F_uv_dx_ini) * T_boost' ; 
F_uv_new_dy_ini = (T_boost * F_uv_dy_ini) * T_boost' ; 
F_uv_new_dz_ini = (T_boost * F_uv_dz_ini) * T_boost' ; 


% Cambiar las coordenadas... más fácil con vectores... con matrices implica
% transformar un vector de matrices... el spinor... Al final hay que hacer
% un boost a diferentes vectores (1,0,0,0) para el tiempo, (0,1,0,0) para
% la X, etc, y en vez de considerar cada vector una unidad, usar las
% matrices F_uv_new_dt_ini, F_uv_new_dx_ini, etc...

[E_dt_aux, B_dt_aux] = SacaCamposTensorEM(F_uv_new_dt_ini);
[E_dx_aux, B_dx_aux] = SacaCamposTensorEM(F_uv_new_dx_ini);

% 13/05/2018.
% Si desde lab avanzo en tiempo dt y en espacio -0.4dx el resultado debería
% ser mismo campo porque es como si el origen no hubiera avanzado nada en
% espacio y en el SRI origen los componentes dtEt son todos cero 
% (el menos es esencial porque lab ve al origen moverse a menos boost!). 
% Luego dtEt – vxdxEx = 0 -> dtEt = vx dxEx (al revés que cómo lo hice)
% Funciona al contrario que como transformaría el tiempo y el espacio
% tengo que repasarme la teoría de boosts en variables covariantes 
% y contravariantes. 
E_dt_test = Sigma * (E_dt_aux + v_boost(1) * E_dx_aux );
B_dt_test = Sigma * (B_dt_aux + v_boost(1) * B_dx_aux );
E_dx_test = Sigma * (E_dx_aux + v_boost(1) * E_dt_aux );
B_dx_test = Sigma * (B_dx_aux + v_boost(1) * B_dt_aux );

% Estos boost indican cómo se combinarán las matrices obtenidas en los ejes
% del SRI origen cuando se pasen a los ejes del SRI destino
Boost_dt = Boost( [1,0,0,0], v_boost);
Boost_dx = Boost( [0,1,0,0], v_boost);
Boost_dy = Boost( [0,0,1,0], v_boost);
Boost_dz = Boost( [0,0,0,1], v_boost);

% Montamos las matrices con los campos en las direcciones t,x,y,z del
% sistema destino
%F_uv_new_dt = F_uv_new_dt_ini * Boost_dt(1) + F_uv_new_dx_ini * Boost_dt(2) + F_uv_new_dy_ini * Boost_dt(3) + F_uv_new_dz_ini * Boost_dt(4);
%F_uv_new_dx = F_uv_new_dt_ini * Boost_dx(1) + F_uv_new_dx_ini * Boost_dx(2) + F_uv_new_dy_ini * Boost_dx(3) + F_uv_new_dz_ini * Boost_dx(4);
%F_uv_new_dy = F_uv_new_dt_ini * Boost_dy(1) + F_uv_new_dx_ini * Boost_dy(2) + F_uv_new_dy_ini * Boost_dy(3) + F_uv_new_dz_ini * Boost_dy(4);
%F_uv_new_dz = F_uv_new_dt_ini * Boost_dz(1) + F_uv_new_dx_ini * Boost_dz(2) + F_uv_new_dy_ini * Boost_dz(3) + F_uv_new_dz_ini * Boost_dz(4);

F_uv_new_dt = F_uv_new_dt_ini * Boost_dt(1) - F_uv_new_dx_ini * Boost_dt(2) - F_uv_new_dy_ini * Boost_dt(3) - F_uv_new_dz_ini * Boost_dt(4);
F_uv_new_dx = - F_uv_new_dt_ini * Boost_dx(1) + F_uv_new_dx_ini * Boost_dx(2) + F_uv_new_dy_ini * Boost_dx(3) + F_uv_new_dz_ini * Boost_dx(4);
F_uv_new_dy = - F_uv_new_dt_ini * Boost_dy(1) + F_uv_new_dx_ini * Boost_dy(2) + F_uv_new_dy_ini * Boost_dy(3) + F_uv_new_dz_ini * Boost_dy(4);
F_uv_new_dz = - F_uv_new_dt_ini * Boost_dz(1) + F_uv_new_dx_ini * Boost_dz(2) + F_uv_new_dy_ini * Boost_dz(3) + F_uv_new_dz_ini * Boost_dz(4);


% Inicializamos las matrices de salida
E_uv_new = zeros(4);
B_uv_new = zeros(4);


% Sacamos los vectores para montar las matrices de derivadas.
[E_dt_new, B_dt_new] = SacaCamposTensorEM(F_uv_new_dt);
[E_dx_new, B_dx_new] = SacaCamposTensorEM(F_uv_new_dx);
[E_dy_new, B_dy_new] = SacaCamposTensorEM(F_uv_new_dy);
[E_dz_new, B_dz_new] = SacaCamposTensorEM(F_uv_new_dz);


E_uv_new = [[0, E_dt_new]; [0,E_dx_new]; [0,E_dy_new]; [0,E_dz_new]];
B_uv_new = [[0, B_dt_new]; [0,B_dx_new]; [0,B_dy_new]; [0,B_dz_new]];



end

