
% Test de la prueba de fuego (crucial)

% Velocidades normalizadas a la velocidad de la luz

% Caso 1 (no empiezo por una dimension porque eso es la excepcionalidad
R =  [1, 1, 0];
v1 = [0, 0, 0];
v2 = [0.4, 0, 0];
% Fallo...

% Caso 2 Si quien se moviese fuese 1...
% R =  [1, 1, 0];
% v1 = [0.4, 0, 0];
% v2 = [0, 0, 0];
% OK


% Caso 3 siguiendo la limintacion anterior
% R =  [1, 1, 2];
% v1 = [0, 0.3, -0.16];
% v2 = [0, 0, 0];
% OK


[Var1, Var2, Var3] = PruebaFuego1(R, v1, v2)