
% NOTA: H es completamente inecesario, por eso siempre valdra 0.
H = [0, 0, 0];

% Juego 1:  OK!!
% ang = pi / 6;
% E_i = [1, 2, 5]


% Juego 2: OK!!
% ang = 0.9 * pi;
% E_i = [-3, 2, 1.4]

 
 % Juego 3: OK!!
% ang = 0.15 * pi;
% E_i = [-3, 0.9, 0]

 
 % Juego 4: OK!!
% ang = 0.25 * pi;
% E_i = [0, 4, -4] 
 

  % Juego 5: OK!!
% ang = 0.5 * pi;
% E_i = [1, 0, 0] 
 
 
 % Juego 6: OK!!
% ang = 0.3 * pi;
% H = [1,2, -3.5];
% E_i = [4, -1.5, 1];

% Juego 7: OK!!
% ang = -0.3 * pi;
%  H = [1.4, -1, 0.9];
% E_i = [2.3, 1.5, -2];

% Juego 8. Da el mismo valor que el inicial
% ang = pi;
% E_i = [1, 0, 0];

% Juego 9. Da el mismo valor que el inicial
% ang = pi;
% E_i = [1, 1, 0];
 
% Juego 10. Cambiando pi/2 invertimos la direccion
% ang = pi/2;
% E_i = [1, 1, 0];

% Juego 11. 
% ang = pi/2;
% E_i = [sqrt(3)/2, 0.5, 0];
 
% Juego 12
% ang = pi/4;
%  E_i = [1, 0, 0];
 
% Juego 13
% ang = pi/8;
% E_i = [1, 0, 0];

% Juego 14. Ver si rota igual con E y H. OK
% E_i = [1, 0, 0];
% H = [0.3, 0.3, 0]
% ang = pi/4;

% Juego 15. OK
%E_i = [1, -2, 0];
%H = [0.5, 0.3, 0]
%ang = pi/2;

% Juego 16. OK
%E_i = [1, -2, 0];
%H = [0.5, 0.3, 0]
%ang = pi;

% Juego 16.
%E_i = [1, -3, 0];
%H = [-0.5, 0.3, 0]
%ang = pi * 20 / 180;

%Juego 17. Con campo en Z
% Giran las componentes Z de Jx y Jy asi como las componentes X e Y de Jz
E_i = [1, -3, 2.3];
ang = pi;

%Juego 18. Con campo en Z
% Giran las componentes Z de Jx y Jy asi como las componentes X e Y de Jz
E_i = [1, -3, 2.3];
H = [0.4, 1.2, 1.7];
ang = pi;

% Valor inicial
Js_i = Corriente_del_momento(E_i, H)

E_rot = Rota_puto_vector(E_i, ang)
H_rot = Rota_puto_vector(H, ang)


Js_o1 = Corriente_del_momento(E_rot, H_rot)

Js_o2 = Rota_corriente_momento(Js_i, ang)