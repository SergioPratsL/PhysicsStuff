
function [e1, e2] = gPerpendicular( V )

dims = size(V);

if dims(1) ~= 1 && dims(2) ~= 3
    disp('mal')
    return
end

x = [1 0 0];
y = [0 1 0];
z = [0 0 1];

v1 = cross(V, x);

if norm(v1) == 0
    v1 = cross(V, y);
end

e1  = v1 / norm(v1);

v2 = cross(V, e1);
e2 = v2 / norm(v2);


v3 = cross( e1, e2);
% No deberia ser necesario pero bueno...
paral = v3 / norm(v3);  

% 21.04.2014. Aplicar regla de la mano derecha
val_aux = dot(paral, V);

% Con comprobar esto debe bastar
if val_aux < 0
    
% Cambiamos el signo de un vector, por lo que he visto deberia bastar 
    e1 = - e1;
end