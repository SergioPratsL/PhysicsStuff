% ¡Aquí comienzan las batallas decisivas de esta campaña!
% Si bien mi idea de este fichero es mostrar como sin el nuevo término, el
% oscilador no puede avanzar.

% 01.01.2026. Esto acabó siendo un cajón de sastre :)

clear;

%k = 0.5; j = [1.2500         0         0         0]
k = 2; 
j = [5     0     0     0];
Phi = [1, 0, 1i*k, 0];


% j  = Obten4CorrienteBispinor(Phi);
% j_norm = j / j(1);
% sigmaX = MatrizPauli(1);
% sigmaX_quad = sigmaX * sigmaX;
% 
% syms x E m w;
% syms y(x);

%eqn = diff(y,x,2) == (E^2 - m^2 + 2*E^2*w^2*x^2 + w^4*x^2)*y

%S = dsolve(eqn)
% (C1*whittakerM(-(E^2 - m^2)/(4*(w^4 + 2*E^2*w^2)^(1/2)), 1/4, x^2*(2*E^2*w^2 + w^4)^(1/2)))/x^(1/2) + (C2*whittakerW(-(E^2 - m^2)/(4*(w^4 + 2*E^2*w^2)^(1/2)), 1/4, x^2*(2*E^2*w^2 + w^4)^(1/2)))/x^(1/2)

%cond = y(0) == 1;
%ySol(x) = dsolve(eqn,cond)
% (C1*whittakerM(-(E^2 - m^2)/(4*(w^4 + 2*E^2*w^2)^(1/2)), 1/4, x^2*(2*E^2*w^2 + w^4)^(1/2)))/x^(1/2) + (gamma((E^2 - m^2)/(4*(w^4 + 2*E^2*w^2)^(1/2)) + 3/4)*whittakerW(-(E^2 - m^2)/(4*(w^4 + 2*E^2*w^2)^(1/2)), 1/4, x^2*(2*E^2*w^2 + w^4)^(1/2)))/(x^(1/2)*pi^(1/2)*(w^4 + 2*E^2*w^2)^(1/8))


[sx, sy, sz] = MatricesPauli();

[gt, gx, gy, gz] = MatricesGamma();
[at, ax, ay, az] = MatricesAlfa();
g5 = MatrizGamma(5);

Sx = MatrizSpin_4_4(1);
Sy = MatrizSpin_4_4(2);
Sz = MatrizSpin_4_4(3);

% Problema, estas matrices no eran, lo del 12.09 está mal:

% matriz_Z = gt * g5 * gy;
% Sy = MatrizSpin_4_4(2);
% % matriz_Z = -Sy
% 
% matriz_Y = gt * g5 * gx
% Sx = MatrizSpin_4_4(1)
% matriz_Y = - Sx
% 
% 
% 
% %C Comprobado:
% %spin_dir = - gt * g5 * g_dir
% matriz_x = -gt * g5 * gx;
% Sx = MatrizSpin_4_4(1);
% 
% matriz_y = - gt * g5 * gy;
% Sy = MatrizSpin_4_4(2);
% 
% matriz_z = - gt * g5 * gz;
% Sz = MatrizSpin_4_4(3);

% LeviCitiva va del 1 al 4
signo_Ay = LeviCivita([3,4,1,2]);
signo_Az = LeviCivita([4,3,1,2]);


% EStas matrices ya tienen el signo.
matriz_y = -1i * gy * g5 * gz;
matriz_z = 1i * gz * g5 * gy;

val_y = ay * matriz_y;
val_z = az * matriz_z;

val_sum = val_y + val_z;

mini_val_z = val_z(1:2,1:2);
mini_val_z_2 = 1i * MatrizPauli(2);

mini_val_y = val_y(1:2,1:2);
mini_val_y_2 = -1i * MatrizPauli(3);


% Esto hace que el Matlab pete!
% syms x E m w k;
% syms y(x);
% 
% eqn =  (E^2 - m^2 + w^4*x^4 -2*E*w^2*x^2 -2*1i*k*(1-2*1i*k*x^2))*y == -diff(y,x,2) + 4*1i*k*x*diff(y,x)
% 
% S = dsolve(eqn)

%syms x w E k m
%syms z(x)

%eqn = diff(z,x,2) == ( w^4*x^4 - 2*E*w^2 + 8*k^2*x^2 + E^2 - m^2)*z
%S = dsolve(eqn)

%x = 1000


spinor_base = [1, 0];   %Z
v = [0, 0.4, 0];
p = fGamma(v) * v;
E = sqrt(1+norm(p)^2);
  
bispinor = DiracSpinorPlainWave(p, spinor_base);
j  = Obten4CorrienteBispinor(bispinor);

bispinor2 = bispinor + [0,0, 0.1i, 0];
j2  = Obten4CorrienteBispinor(bispinor2);
bispinor3 = bispinor + [0,0, 0.1, 0];
j3  = Obten4CorrienteBispinor(bispinor3);


% 01.01.2026
signo_Ay_dx = LeviCivita([1,4,3,2]);
nueva_matriz = g5 * gt * gz * signo_Ay_dx;
producto_spinor_Z = nueva_matriz * [1,0,0,0]';

% 25.01.2026
signo_V_dx = LeviCivita([3,4,2,1]);
matriz_Ex = g5 * gy * gz * signo_V_dx


% Ecuaciones diferenciales (fueron una mierda)

%syms y(t) z(t)
%eqns = [diff(y,t) == z, diff(z,t) == -y];
%S = dsolve(eqns)


%syms E w k x
%syms f(x) g(x)
% Con fallos
%eqns = [diff(f,x) == (E - w*x  - 1i*k*x)*f - k*x*g, diff(g,x) == -(E - w*x + 1i*k*x)*g - k*x*f];
% Completa
%eqns = [diff(g,x) == (E + w*x^2  - 1i*k*x)*f - k*x*g, diff(f,x) == -(E + w*x^2 + 1i*k*x)*g + k*x*f];

% Sin constantes
%eqns = [diff(g,x) == (1 + x^2  - 0.001*1i*x)*f - 0.001*x*g, diff(f,x) == -(1 + x^2 + 0.001*1i*x)*g + 0.001*x*f];

% Esta se resolvió!
%eqns = diff(f,x) == (1 + x^2  - 1*1i*x)*f;

% Sin nuevos términos
%eqns = [diff(g,x) == (E + w*x^2)*f, diff(f,x) == -(E + w*x^2)*g];

%S = dsolve(eqns)


%syms E a B x
%syms f(x) g(x)
%eqns = [diff(g,x) == (E/a)*f- a*B*x*g, diff(f,x) == (E/a)*g+a*B*x*g];

%S = dsolve(eqns)


% Una mierda como de costumbre, puto Matlab.
%syms E m B x
%syms f(x) g(x)
%eqns = [diff(f,x) == B*x*f + B*g, diff(g,x) == -(1/B)*(E-m)^2*f + B*x*g];
%S = dsolve(eqns)