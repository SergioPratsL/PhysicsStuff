clear;

% Hay que modificar el GetThomasMatrix para quitar el ; que muestra el
% ángulo o no hay manera :(.

vBA = [0.4, 0, 0];
%vBA_new = [0.4, 0.0001, 0];    % Perfecto
%vBA_new = [0.4002, 0.0001, 0]; % Casi perfecto ... ratio=0.9997
%vBA_new = [0.40002, 0.0001, 0]; % ratio de 0.998
%vBA_new = [0.40002, 0.001, 0]; % ratio de 1.0000
vBA_new = [0.400002, 0.0001, 0]; % ratio de 1.0000

% Visto desde el lab!
dv = vBA_new - vBA;

[rotMatrixGetThomasRotMatrix, anguloGetThomasMatrix] = GetThomasRotMatrix(vBA, vBA_new);
anguloGetThomasMatrix = anguloGetThomasMatrix

[eje_rot, angulo] = ThomasPrecessionAng(0, dv, vBA);
anguloThomasPrecession = angulo

ratio = anguloGetThomasMatrix / anguloThomasPrecession





