function F = FLorentz(v, E, B)

% Fuerza de Lorentz sin ninguna unidad
% velocidad normalizada
% carga = 1

F = E + cross(v, B);

