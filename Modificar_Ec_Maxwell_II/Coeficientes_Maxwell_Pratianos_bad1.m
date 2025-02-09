function [ dt_E, dt_B ] = Coeficientes_Maxwell_Pratianos_bad1( E, B, a, b )
% Esta funci�n debe calcular la variaci�n del campo electromagn�tico
% necesaria para que con dicha variaci�n se pueda compensar la variaci�n 
% de energ�a y momento del campo radiado aisalo de la(s) part�culas

% Las variables E y B representan el campo el�ctrico TOTAL que hay en un
% punto. 
% La variable a es la ganancia o p�rdida de energ�a del campo radiado de 
% la part�cula estudiada (es decir, su campo radiado aislado). Por ganancia
% o p�rdida de energ�a se entiende no la evoluci�n temporal de la energ�a
% sino dt_u  + div(S) = -a
% An�logamente b es la ganancia o p�rdida de momento del campo radiado 
% de la part�cula estudiada y se entiende como dt_S - Div(TensMaxwell) = b
% PD: el por qu� del -Div() en vez +Div de no lo s�

% Hay que resolver un sistema de ecuaciones cuyas incognitas son 
% [Ey, Ez, By, Bz], despues se obtienen Ex y Bx como
% Ex = - Ey - Ez, Bx = -By - Bz
% Esto se deriva de aplicar la ecuaci�n de Gauss a la variaci�n del campo
% debida al EM y a que el efecto de la carga ya aplic� al campo "normal".


% Preparar los coeficientes del sistema de ecuaciones
Cofs = zeros(4);

Valores = [(-2*a) (-b)]';

Cofs(1,1) = E(2) - E(1);
Cofs(1,2) = E(3) - E(1);
Cofs(1,3) = B(2) - B(1);
Cofs(1,4) = B(3) - B(1);

Cofs(2,1) = B(3);
Cofs(2,2) = -B(2);
Cofs(2,3) = -E(3);
Cofs(2,4) = E(2);

Cofs(3,1) = B(3);
Cofs(3,2) = B(1) + B(3);
Cofs(3,3) = -E(3);
Cofs(3,4) = -E(1) - E(3);

Cofs(4,1) = -B(1) - B(2);
Cofs(4,2) = -B(2);
Cofs(4,3) = E(1) + E(2);
Cofs(4,4) = E(2);


Result = Cofs \ Valores;

dt_E = [- Result(1) - Result(2), Result(1), Result(2)];
dt_B = [- Result(3) - Result(4), Result(3), Result(4)];

end

