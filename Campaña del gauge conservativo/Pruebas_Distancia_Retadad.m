
% Prueba 1: todos sobre el eje X. OK
%R = [2, 0, 0];
%v = [0.4, 0, 0];

% Prueba 2. OK
% R = [-2, 0, 0];
% v = [0.4, 0, 0];

% Prueba 3. OK
% R = [2, 0, -1];
% v = [0.25, -0.2, 0];

% Prueba 4. OK
% R = [0, 1, -1.5];
% v = [0, -0.2, 0.3];

% Prueba 5. OK
R = [10, -14, 6];
v = [0.29, -0.24, 0.18];
% Rr = 18.9373  -21.3964   11.5473

% Parece que todas las pruebas dan bien, realmente me he basado en que
% satisfagan las 4 ecuaciones de la entrada del 06.10.2015 del diario de
% fisica

Rr = DRetardada(R, v);


% Aprovecho el script...

% Prueba 5. OK
R = [1, -2, 0];
v = [0.29, -0.24, 0.18];

R_ret = PostRetFut(R, 0.3, 'r')
% Hay que acerlo asi para que vaya bienm si inviertes en el paso anterior
% la lias y obtienes el vector futuro pero invertido.
R_ret = - R_ret;

R_fut = PostRetFut(R, 0.3, 'f')



