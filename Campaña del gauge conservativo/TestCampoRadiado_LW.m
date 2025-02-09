
% Prueba 1. OK
Rx = 0;
Ry = 1;
a = [1, 0, 0];


% Prueba 2. OK
%Rx = 0;
%Ry = 1.4;
%a = [0, 0, 1];

% Prueba 3. OK
%Rx = -1.1;
%Ry = 1.4;
%a = [0, 1, 1];


% Prueba 4. OK
% Rx = -1.9;
% Ry = 1.4;
% a = [2, 0, 1];

% Prueba 5. OK
% Rx = -1.9;
% Ry = 0;
% a = [-1, 1, 1];

 % Prueba 6. OK
% Rx = -1.9;
% Ry = 0.76;
% a = [0, 0, -2.1];


% Prueba 7. OK
 Rx = 1;
 Ry = 1;
 a = [-1, -1, 0];

 % Prueba 8. OK
 Rx = 1;
 Ry = 1;
 a = [-1, -1, 3];
 
 
 
v_0 = [0, 0, 0];

R = [Rx, Ry, 0];

[E_deriv_LW, B_deriv_LW] = CampoRadiadoDerivadas_LW( Rx, Ry, a)

[E, B] = CampoRadiadoCargaAcelerada(R, v_0, a)

Dif_E = E - E_deriv_LW 

Dif_B = B - B_deriv_LW


