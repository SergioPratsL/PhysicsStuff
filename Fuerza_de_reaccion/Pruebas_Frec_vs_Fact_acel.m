% El objetivo es ver c�mo queda la componente avanzada del campo radiado
clear;      

% Mi unidad de tiempo ser� el microsegundo.
t_unit = 10^-6;
c = 299792458 * t_unit;
% Mi unidad de distancia ser� el tiempo que la luz en un microsegundo
d_unit = c;

% distancia base representa un metro en las unidades de d_unit...
r_base = 1 / d_unit;

% La aceleraci�n se hace en base a t_unit.

% r representar� la posici�n actual del objeto, usada para determinar las
% posiciones pasadas y futura.


% Prueba 1 - Campos inducidos E y B se compensan, campo radiado E tb, B no.
% velocidad a la que se mueve el observador respecto a la carga 
%v = [0, 0.2, 0];
v = [0, 0, 0.4];
%v = [0, 0, 0];
% posici�n relativa a la carga del punto en el que evaluamos el campo
r = [1, 0, 0];   % esto ser�an r_base metros.
%r = r_base * [1, 0, 0];   
% aceleracion propia
%a = 0.4 * [0, 0, 1];
a = 10^-1 * [0, 0, 1];
%a = 10^-3 * [0, 0, 1];


% Nos centramos en el receptor, calcular la posici�n y velocidad pasadas 
% y futura, como primera aproximaci�n para la posici�n podemos asumir 
% velocidad constante.

% https://en.wikipedia.org/wiki/Retarded_time
[rt_ret, v_ret] = Obten_Posicion_Ret_O_Avz(r, v, a, 'r');
[rt_avz, v_avz] = Obten_Posicion_Ret_O_Avz(r, v, a, 'a');

[Ei_ret, Bi_ret, Er_ret, Br_ret] = Campo_EM_Total_No_Units(rt_ret, v_ret, a);

% La aceleraci�n no cambia de signo
[Ei_avz, Bi_avz, Er_avz, Br_avz] = Campo_EM_Total_No_Units(rt_avz, -v_avz, a);

Ei_ret
Ei_avz;


E_ret = Ei_ret + Er_ret;
B_ret = Bi_ret + Br_ret;
E_Dif_tot = Ei_ret + Er_ret - Ei_avz - Er_avz;
B_Dif_tot = Bi_ret + Br_ret + Bi_avz + Br_avz;


Fret = FLorentz(v, E_ret, B_ret)
Ftot = FLorentz(v, E_Dif_tot, B_Dif_tot)

% Veamos la fuerza que a�n se hace, proporcionalmente.
ratio = norm(Ftot) / norm(Fret)
