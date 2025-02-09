function [V, A] = PotencialDeYinnAmpliado(R, v)

% La diferencia entre este potencial y el script "PotencialDeYinn" es que
% aqui la velocidad puede ser en cualquier direccion

% Es el potencial desde el SRI del emisor,
% v es la velocidad del receptor (vectorial)
% R es el vector retardado que va del emisor al receptor (sin tiempo)

if norm(v) == 0
    return
end 

if norm(R) == 0        % :@
    return
end

distance = norm(R);
speed = norm(v);

V = 1 / distance;

v_x_r = cross(v, R);

r_x_v_r = cross( R, v_x_r);
v_x_v_r = cross( v, v_x_r );

X1 = - r_x_v_r  / dot(R, v_x_v_r) / distance;

X2 = v_x_v_r / dot(v_x_r, v_x_r) / speed;

A = X1 + X2;

