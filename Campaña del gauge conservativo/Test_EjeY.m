
% Prueba piloto, esta se supera finalmente
% OK a la primera (ni un error de compilacion, me quedo estupefacto) 
% --> Los cojones a la primera
% Rx = 0;
% Ry = 1;
% v = 0.4;

% Otra prueba. OK
%  Rx = 1;
%  Ry = 1;
%  v = 0.3;


% Prueba 3: OK
% Rx = -2;
% Ry = 1.4;
% v = 0.35;

% Prueba 4: OK
% Rx = 1;
% Ry = 3;
% v = -0.6;

% Prueba 5: OK
% Rx = 3.3;
% Ry = -2;
% v = 0.99;

% Prueba 6: OK
% Rx = 0.6;
% Ry = -1.1;
% v = 0.01;

% Prueba 7: 
%  Rx = -100;
%  Ry = 1;
%  v = 0.4;

 
 % Prueba para campo radiado
Rx = 0;
Ry = 1;
v = 0.4;
 

[Pot_A, Pot_B, F_rad, Dif] = DifPot_EjeY( Rx, Ry,  v)