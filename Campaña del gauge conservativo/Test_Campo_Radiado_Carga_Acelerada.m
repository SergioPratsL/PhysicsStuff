
% Prueba 1 (OK):
% R = [10, 5, 4];
% v = [0, 0, 0];
% dv = [0.2, 0.1, 0.08];
% Resultado:
% Eo = [     0     0     0 ]
% Bo = [     0     0     0 ]


% Prueba 2 (OK):
% R = [10, 5, 4];
% v = [0, 0, 0];
% dv = [-0.1, 0.1, 0.125];
% Resultado:
% Eo = [ 0.0084   -0.0084   -0.0105 ]
% Ho = [ -0.0016    0.0117   -0.0106 ]
% Varias verificacoines hechas: Eo paral a dv (perpend a R), Eo es dv / R.
% B es igual que E pero en direccion ortogonal, y tambien ortogonal a V


% Prueba 3 (OK):
% R = [1, 0, 0];
% v = [0, 0.4, 0];
% dv = [0, 0, 1];
% Resultados:
% Eo = [ 0         0   -0.8004 ]
% Bo = [ 0    0.8004         0 ]
% La reducion respecto a 1 es debida al aumento de la distancia actual.

% Tiene buena pinta, aunque no pondria la mano en el fuego.

clear

R = [1, 0, 0];
v = [0, 0.1, 0];
dv = [0, 0, 1];

[Eo, Bo] = CampoRadiadoCargaAcelerada( R, v, dv)