% Matlab para verificar que en partículas inerciales el campo futuro y el
% pasado han de ser iguales visto desde cualquier SRI.
% Es obvio que en el SRI de la carga Fret-Fact = 0, la idea es confirmar
% que en otros SRIs también es asço

clear;      

% Prueba 1 - OK
% velocidad a la que se mueve el observador respecto a la carga
%v = [0.4, 0, 0];
% posición relativa a la carga del punto en el que evaluamos el campo
%r = [1, 0, 0];

% Prueba 2 - Campo eléctrico bien, campo magnético tiene signo opuesto!!!
v = [0.4, 0, 0];
r = [0, 1, 0];


% Prueba 3 - Campo eléctrico bien, campo magnético tiene signo opuesto!!! 
%v = [0.4, 0, 0];
%r = [-1.4, 1, 0];


% Prueba 4 - Campo eléctrico bien, campo magnético tiene signo opuesto!!!
%v = [0.2, 0, -0.34];
%r = [-0.4, 1.2, 0.6];


rt_ret = [norm(r), r];
rt_avz = [-norm(r), r];

% lo que acaba en 2 está visto por el otro SRI.
rt2_ret = Boost(rt_ret, v);
r2_ret = rt2_ret(2:4)

rt2_avz = Boost(rt_avz, v);
r2_avz = rt2_avz(2:4)

[E2_ret, B2_ret] = CampoInducido_sin_unidades(r2_ret, -v)

[E2_avz, B2_avz] = CampoInducido_sin_unidades(r2_avz, v)

Dif_E = E2_ret - E2_avz

% Por lo visto se ha de invertir el campo magnético porque la evolución
% hacia atrás en el tiempo requiere invertir el signo, en el eléctrico no
% pasa ya que se invierte sobre el magnético que ya está invertido.
Dif_B = B2_ret - (- B2_avz)




