
clear 

% La madre de todas las batallas!!

% Caso 1 auxiliar
Rx = 0;
Ry = 1;
v = 0.4;
t0_B = 0;
a_A = [-1,0,0];
a_B = [1,0,0];


% Caso de mierda
%Rx = 1;
%Ry = 1;
%v = 0.3;
%t0_B = 0;
%a_A = [-1,0,0];
%a_B = [1,0,0];


% Caso de mierda
%Rx = 2;
%Ry = -1;
%v = 0.3;
%t0_B = 0;
%a_A = [-1,0,0];
%a_B = [1,0,0];




[Dif_Pot_A, Dif_Pot_B, Dif_Pot_Tot, Dif_Pot_Tot_en_B] = DifPotAcelDual(Rx, Ry, v, t0_B, a_A, a_B)
