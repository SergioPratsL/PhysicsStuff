clear all; clc;

% --- Physical Constants (SI Units) ---
m = 9.10938356e-31;     % Electron mass (kg)
c = 299792458;          % Speed of light (m/s)
hbar = 1.0545718e-34;   % Reduced Planck constant (J.s)
q = 1.602176 * 10^-19;


spinor_base = [1, 0];   %Z
% USANDO DiracSpinorPlainWave_GtGt
% Con v=.01 el spin predominante es el segundo como era de esperar, pero si
% la velocidad es en X o Y parece que la componente es la 4ª, mientras que
% si es en Z la componente es la 3, es decir, para px y py el spin se
% invierte, para pz no.
% ~USANDO DiracSpinorPlainWave_GtGt

% Usando la versión mejorada DiracSpinorPlainWave_GtGtSp1,
% los spinors son idénticos que con el Hamiltoniano oficial, 
% eso no era necesario, bastaba que las corrientes fueran la
% correctas, pero en todo caso ¡PRUEBA SUPERADA!

%v = [0.01, 0, 0];
%v = [0, 0.01, 0];
%v = [0, 0, 0.01];

% Las corrientes están bien, pero "Phi" no es más que un pelele.
%v = [0.4, 0, 0];
%v = [0, 0.4, 0];
%v = [0, 0, 0.4];

%v = [0, 0, 0];

v = [0.5, -0.37, 0.4];

p = fGamma(v) * v;
E = sqrt(1+norm(p)^2);
  
% Legado
%bispinor = DiracSpinorPlainWave_GtGt(p, spinor_base).'

% Este usa las fórmuls que se basan en el Hamiltoniano oficial, no usan el
% gt para la masa.
bispinor_normal = DiracSpinorPlainWave(p, spinor_base);
bispinor_normal = bispinor_normal / norm(bispinor_normal)

% Este se basa en mi modificación del Hamilt.
bispinor_2 = DiracSpinorPlainWave_GtGtSp1(p, spinor_base)



[jt, jx, jy, jz]  = ObtenCorrientesBispinor(bispinor_2);
corrientes_intrinsecas = [jt, jx, jy, jz] / jt