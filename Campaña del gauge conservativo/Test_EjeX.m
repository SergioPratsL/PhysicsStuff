
% 09.06.2016. Aunque este script es del 12.2016, pruebo ahora para
% compararlo con la version adoptada a la nueva formulacion, que falla...

% Prueba piloto, esta se supera finalmente
% Rx = 0;
% Ry = 1;
% v = 0.4;

% Otra prueba, a ver si se observa un "patron de error
% Rx = 1;
% Ry = 1;
% v = 0.3;
% Tb OK

% Prueba 3: OK!!!
% Rx = -2;
% Ry = 1;
% v = 0.3;

% Prueba 4: OK!!!
Rx = 1;
Ry = 2;
v = -0.6;

% Prueba 5: OK!!
% Rx = 3.3;
% Ry = -2;
% v = 0.99;

% Prueba 6: 
%  Rx = 1.4;
%  Ry = -1.1;
%  v = 0.01;

%OK, KOKOKOKOKOKOKOKOKOK!!!!! Otra victoria :D:D:D:D:D:D:D:D
 

% Prueba para comparar con campo radiado
% Prueba 1.
% Rx = 0;
% Ry = 1;
% v = 0.3;
%a = [1, 0, 0];


[Pot_A, Pot_B, F_rad, Dif] = DifPot_EjeX( Rx, Ry,  v)