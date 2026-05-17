
% ESTO ESTA MAL

% (E-1)x=My
% (E-1)y=Mx
% (E-1)x=M*Mx/(E-1)
% (E-1)^2*x=M^2*x
% Pero esto no es verdad, porque el M^2 se vuelve una identidad, luego es motivo de cese.

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
% correctas, pero en todo caso 
% PERO EN TODO CASO ERA UNA FARSA, HICE UNA PIRULA QUE LO INVALIDA TODO.

%v = [0.01, 0, 0];
%v = [0, 0.01, 0];
%v = [0, 0, 0.01];

% Las corrientes están bien, pero "Phi" no es más que un pelele.
%v = [0.4, 0, 0];
%v = [0, 0.4, 0];
v = [0, 0, 0.4];

%v = [0, 0, 0];

%v = [0.5, -0.37, 0.4];

p = fGamma(v) * v;
E = sqrt(1+norm(p)^2);
  
% Legado
%bispinor = DiracSpinorPlainWave_GtGt(p, spinor_base).'

% Este usa las fórmuls que se basan en el Hamiltoniano oficial, no usan el
% gt para la masa.
bispinor_normal = DiracSpinorPlainWave(p, spinor_base);
bispinor_normal = bispinor_normal / norm(bispinor_normal)

% Este se basa en mi modificación del Hamilt.
%bispinor_2 = DiracSpinorPlainWave_GtGtSp1(p, spinor_base)

bispinor_2 = DiracSpinorPlainWave_GtGt(p, spinor_base)


[jt, jx, jy, jz]  = ObtenCorrientesBispinor(bispinor_2);
corrientes_intrinsecas = [jt, jx, jy, jz] / jt;


E = sqrt(norm(p)^2 + 1);
P_Dot_PauliVector = PauliVectorEscalarProd(p);


% Con el Hamiltoniano normal, H_bispinor = H_esperado
% Phi = bispinor_normal(1:2);
% Chi = bispinor_normal(3:4);
% 
% H_Phi = Phi + (P_Dot_PauliVector * Chi.').';
% % Aquí el -Chi está bien porque la ecuación tiene -m
% H_Chi = -Chi + (P_Dot_PauliVector * Phi.').';
% 
% H_bispinor = [H_Phi, H_Chi]
% H_esperado = E * bispinor_normal


% Phi = bispinor_2(1:2);
% Chi = bispinor_2(3:4);
% 
% P_Dot_PauliVector = PauliVectorEscalarProd(p);
% 
% 
% H_Phi = Phi + (P_Dot_PauliVector * Chi.').';
% % Hack
% %H_Phi = - Phi + (P_Dot_PauliVector * Chi.').';
% H_Chi = Chi + (P_Dot_PauliVector * Phi.').';
% 
% H_bispinor = [H_Phi, H_Chi]
% 
% 
% H_esperado = E * bispinor_2
% 
% 
% %  ratio_1 = SmartRatio(H_bispinor(1), bispinor_2(1))
% %  ratio_2 = SmartRatio(H_bispinor(2), bispinor_2(2))
% %  ratio_3 = SmartRatio(H_bispinor(3), bispinor_2(3))
% %  ratio_4 = SmartRatio(H_bispinor(4), bispinor_2(4))









