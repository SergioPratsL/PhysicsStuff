% El objetivo es ver cómo queda la componente avanzada del campo radiado
clear;      

% Mi unidad de tiempo será el microsegundo.
t_unit = 10^-6;
c = 299792458 * t_unit;
% Mi unidad de distancia será el tiempo que la luz en un microsegundo
d_unit = c;

% distancia base representa un metro en las unidades de d_unit...
r_base = 1 / d_unit;

% La aceleración se hace en base a t_unit.

% r representará la posición actual del objeto, usada para determinar las
% posiciones pasadas y futura.


% Prueba 1 - Campos inducidos E y B se compensan, campo radiado E tb, B no.
% velocidad a la que se mueve el observador respecto a la carga 
v = [0, 0.2, 0];
% posición relativa a la carga del punto en el que evaluamos el campo
r = r_base * [1, 0, 0];   % esto serían 1 metros.
%r = [1, 0, 0]; 
% aceleracion propia
a = 10^-3 * [0, 0, 1];
%a = 10^-1 * [0, 0, 1];




% Prueba 2 - Aquí los campos radiados son cero así que todo ok
% velocidad a la que se mueve el observador respecto a la carga 
%v = [0, 0.4, 0];
% posición relativa a la carga del punto en el que evaluamos el campo
%r = r_base * [1, 0, 0];
% aceleracion vista desde el laboratorio, no aceleración propia
%a = 10^-3 * [1, 0, 0];


% Prueba 3 - Los campos inducidos se anulan pero los radiados no, ni E ni B
% velocidad a la que se mueve el observador respecto a la carga 
%v = [0, 0.4, 0];
% posición relativa a la carga del punto en el que evaluamos el campo
%r = r_base * [1, 0, 0];
% aceleracion vista desde el laboratorio, no aceleración propia
%a = 10^-3 * [0, 1, 0];


% Prueba 3B - Confrontaremos los campos con igualdad de condiciones.
%v = [0, 0.4, 0];
% posición relativa a la carga del punto en el que evaluamos el campo
%r = r_base * [1, 0, 0];
% aceleracion vista desde el laboratorio, no aceleración propia
%a = 10^-3 * [0, 1, 0];


% Prueba 4B - Campo Bi es residual y sus diferencias son relativamente
% grandes ¿Errores de precisión? campo Ei error orden 10^-5...
% diferencias en campo radiado grandes tanto en B como en E... kk
% velocidad a la que se mueve el observador respecto a la carga 
%v = [0.4, 0, 0];
% posición relativa a la carga del punto en el que evaluamos el campo
%r = r_base * [1, 0, 0];
% aceleracion vista desde el laboratorio, no aceleración propia
%a = 10^-3 * [0, 1, 0];


% Nos centramos en el receptor, calcular la posición y velocidad pasadas 
% y futura, como primera aproximación para la posición podemos asumir 
% velocidad constante.

% https://en.wikipedia.org/wiki/Retarded_time
rt_ret = Posicion_Retardada_Bruta(r, v, a);
rt_avz = Posicion_Avanzada_Bruta(r, v, a);

v_ret = v - a * norm(rt_ret);
v_avz = v + a * norm(rt_avz);

[Ei_ret, Bi_ret, Er_ret, Br_ret] = Campo_EM_Total_No_Units(rt_ret, v_ret, a);

% La aceleración no cambia de signo
[Ei_avz, Bi_avz, Er_avz, Br_avz] = Campo_EM_Total_No_Units(rt_avz, -v_avz, a);


E_ret = Ei_ret + Er_ret;
B_ret = Bi_ret + Br_ret;
E_Dif_tot = Ei_ret + Er_ret - Ei_avz - Er_avz;
B_Dif_tot = Bi_ret + Br_ret + Bi_avz + Br_avz;

Ei_ret


