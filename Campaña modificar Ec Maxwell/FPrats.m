
function F = FPrats(v, E, B)
% Esta es la fuerza que se obtiene de aplicar a la ecuación de Ampere el
% termino pratiano de Ev^2 * J donde Ev es el campo en la direccion de la
% velocidad.

% Aqui se obvian todas las constantes, esta prueba es solo para ver si esta
% fuerza realmente es coherente ante un boost.

% En el caso viciado de v = 0, devolver un cero.
if norm(v) == 0
    F = [0, 0, 0];
    return
end         

v_mod = norm(v);

v_norm = v / v_mod;

% Obtener la maldita base!

[e1, e2, e3] = ObtenBase( v_norm);

%T = CambioBaseEnh(E, v);
%T_inv = T^(-1);

%E_rot = T * E';
%B_rot = T * B';

E_rot = CambioBase(e1, e2, e3, E);
B_rot = CambioBase(e1, e2, e3, B);

F1 = E_rot(1)^3;

F2 = 3 * E_rot(1)^2 * (E_rot(2) - v_mod * B_rot(3));
F3 = 3 * E_rot(1)^2 * (E_rot(3) + v_mod * B_rot(2));
% Cambiar signo (fallo) y hace que pruebas superadas fallen.
%F2 = 3 * E_rot(1)^2 * (E_rot(2) + v_mod * B_rot(3));
%F3 = 3 * E_rot(1)^2 * (E_rot(3) - v_mod * B_rot(2));
% Tambien falla esta, pero las pruebas anteriores se pasan:
%F2 = - 3 * E_rot(1)^2 * (E_rot(2) - v_mod * B_rot(3));
%F3 = - 3 * E_rot(1)^2 * (E_rot(3) + v_mod * B_rot(2));

F_rot = [F1, F2, F3];

F = F1 * e1 + F2 * e2 + F3 * e3;
%F = (T_inv * F_rot')';

