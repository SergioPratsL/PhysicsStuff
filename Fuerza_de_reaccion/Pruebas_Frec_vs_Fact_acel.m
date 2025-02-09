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
%v = [0, 0.2, 0];
v = [0, 0, 0.4];
%v = [0, 0, 0];
% posición relativa a la carga del punto en el que evaluamos el campo
r = [1, 0, 0];   % esto serían r_base metros.
%r = r_base * [1, 0, 0];   
% aceleracion propia
%a = 0.4 * [0, 0, 1];
a = 10^-1 * [0, 0, 1];
%a = 10^-3 * [0, 0, 1];


% Nos centramos en el receptor, calcular la posición y velocidad pasadas 
% y futura, como primera aproximación para la posición podemos asumir 
% velocidad constante.

% https://en.wikipedia.org/wiki/Retarded_time
[rt_ret, v_ret] = Obten_Posicion_Ret_O_Avz(r, v, a, 'r');
[rt_avz, v_avz] = Obten_Posicion_Ret_O_Avz(r, v, a, 'a');

[Ei_ret, Bi_ret, Er_ret, Br_ret] = Campo_EM_Total_No_Units(rt_ret, v_ret, a);

% La aceleración no cambia de signo
[Ei_avz, Bi_avz, Er_avz, Br_avz] = Campo_EM_Total_No_Units(rt_avz, -v_avz, a);

Ei_ret
Ei_avz;


E_ret = Ei_ret + Er_ret;
B_ret = Bi_ret + Br_ret;
E_Dif_tot = Ei_ret + Er_ret - Ei_avz - Er_avz;
B_Dif_tot = Bi_ret + Br_ret + Bi_avz + Br_avz;


Fret = FLorentz(v, E_ret, B_ret)
Ftot = FLorentz(v, E_Dif_tot, B_Dif_tot)

% Veamos la fuerza que aún se hace, proporcionalmente.
ratio = norm(Ftot) / norm(Fret)
