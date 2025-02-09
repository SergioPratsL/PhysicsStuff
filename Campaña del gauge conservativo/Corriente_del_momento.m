function Js = Corriente_del_momento(E, H)

% NOTA: esta funcion no tiene unidades para E y H
% suelo tratar (burdamente) B y H como si fueran lo mismo


Js = zeros(3,3);

E_gir = E';
H_gir = H';

Num = 0.5 * (norm(E)^2 + norm(H)^2);

Js = kron(E, E_gir) + kron(H, H_gir) - eye(3) * Num;

% Formula mala usada en algunas pruebas
%Js = kron(E, E_gir) + kron(H, H_gir);