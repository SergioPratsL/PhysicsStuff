% Script para ver cómo cambia el spin con la rotación, para ver qué es eso
% de que el spin ha de cambiar 720º para volver a su posición inicial.
% 06.01.2026. Este fichero se basa en C:\Fisica\Scripts matlab\HacerTuTeoria\SpinOMuerte\PruebasRotacionSpin.m

% p=0, sistema en reposo!!!
 p0 = [0, 0, 0];


% Rotacion en el eje del spin, hace falta rotar 2pi para volver al valor inicial. 
 
% Sz = MatrizSpin_4_4(3);
% 
% ph_ini = [1, 0];     %  Z
% Psi_ini = DiracSpinorPlainWave(p0, ph_ini)
% Psi = Psi_ini;
% 
% dfase = pi / 1000;
% n = 0;
% %while n < 2000      % [1, 0, 0, 0];  Fase 2pi
% %while n < 1500      % [1i, 0, 0, 0]; Fase 3pi/2
% while n < 1000     % [-1, 0, 0, 0];  Fase pi
% %while n < 750      % [-0.7071 - 0.7071i, 0        0         0]; Fase 3pi/2
% %while n < 500     % [-1i, 0, 0, 0];  Fase pi/2
% %while n < 250     % [0.7071 - 0.7071i, 0        0         0];  Fase pi/4
% 
%     dPsi = dfase * SpinorRotation(Sz, Psi)';
%     Psi = Psi + dPsi;
%     
%     Psi = Psi / norm(Psi);
%     
%     n = n + 1;
% end
% 
% fase_fin = dfase * n
% Psi = Psi


% Otra prueba: ver cómo cambia la fase global rotando en el eje del spin,
% se tarda la mitad de fase (pi) en volver a la misma posición espacial, 
% y lo mismo (2pi) en volver al mismo valor

Sx = MatrizSpin_4_4(1); % X

ph_ini = [1, 0];     %  Z
Psi_ini = DiracSpinorPlainWave(p0, ph_ini)
Psi = Psi_ini;

dfase = pi / 1000;
n = 0;
% while n < 2000     % [1, 0, 0, 0];  Fase 2pi
% %while n < 1500     % [0, 1i, 0, 0];  Fase 3pi/2
% %while n < 1000     % [-1, 0, 0, 0];  Fase pi
% %while n < 500     % [0, -1i, 0, 0];  Fase pi/2
% 
%     dPsi = dfase * SpinorRotation(Sx, Psi)';
%     Psi = Psi + dPsi;
%     
%     Psi = Psi / norm(Psi);
%     
%     n = n + 1;
% end
% 
% fase_fin = dfase * n
% Psi = Psi


dPsi_test = SpinorRotation(Sx, Psi_ini)'
