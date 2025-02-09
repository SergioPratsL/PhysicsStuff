function [ dt_E, dt_B ] = Coeficientes_Maxwell_Pratianos_bad1( E, B, a, b )
% Esta función debe calcular la variación del campo electromagnético
% necesaria para que con dicha variación se pueda compensar la variación 
% de energía y momento del campo radiado aisalo de la(s) partículas

% Las variables E y B representan el campo eléctrico TOTAL que hay en un
% punto. 
% La variable a es la ganancia o pérdida de energía del campo radiado de 
% la partícula estudiada (es decir, su campo radiado aislado). Por ganancia
% o pérdida de energía se entiende no la evolución temporal de la energía
% sino dt_u  + div(S) = -a
% Análogamente b es la ganancia o pérdida de momento del campo radiado 
% de la partícula estudiada y se entiende como dt_S - Div(TensMaxwell) = b
% PD: el por qué del -Div() en vez +Div de no lo sé

% Hay que resolver un sistema de ecuaciones cuyas incognitas son 
% [Ey, Ez, By, Bz], despues se obtienen Ex y Bx como
% Ex = - Ey - Ez, Bx = -By - Bz
% Esto se deriva de aplicar la ecuación de Gauss a la variación del campo
% debida al EM y a que el efecto de la carga ya aplicó al campo "normal".


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

