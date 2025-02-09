function [V, A] = PotencialDeYinn(R, v)

% Potencial obtenido por integracion directa, definitivo!!!

% ...Incluyo la constante que es lo que me parece menos definitivo...

% Es el potencial desde el SRI del emisor,
% v es la velocidad del receptor
% R es la distancia retardada, del emisor al receptor (sin tiempo)

if norm(v) == 0
    return
end 

if norm(R) == 0        % :@
    return
end
 
% Esto es lo que en mis papeles es la X 
speed = norm(v);
v_norm = v / speed;
 
distance = norm(R);
R_norm = R / distance;

e_2_not_norm = R - v_norm * dot(R, v_norm);

vector_ortogonal = cross( v_norm, R_norm );
 
norm_v_ort = norm(vector_ortogonal);
 
if  norm_v_ort > 0
     
    vector_ortonormal = vector_ortogonal / norm(vector_ortogonal);

% En este caso la distancia y la diferencia de lecturas estaban alineadas!!
% elegir el vector que me salga del rabo!
else
    
    if dif_v_norm(2) ~= 0  ||  dif_v_norm(3) ~= 0
        vect_test = [1, 0, 0];
    else
        vect_test = [0, 1, 0];
    end
    
    vect_ortogonal = cross( dif_v_norm, vect_test );
    
    vector_ortonormal = vector_ortogonal / norm(vector_ortogonal);    
    
end


if norm(e_2_not_norm) == 0
    e_2 = cross( v_norm, vector_ortonormal );
    e_2 = e_2 / norm(e_2);      % Por si las moscas
else
    e_2 = e_2_not_norm / norm(e_2_not_norm);
end
    
    
% Ahora puedo trabajar con mis formulas

Rx = dot(R, v_norm);
Ry = dot(R, e_2);

% Potencial en la base calculada:
V = 1 / norm(R);

A_new_base(1) = 1 / (distance * speed );
% A_new_base(1) = A_new_base(1) + 800; % Si es sensible


% Ojo, ojito.
if Ry == 0
    A_new_base(2) = inf;        % :P
    
else   
    
    Constante = - 1 / ( abs(speed)  * Ry );
%    Constante = 9999; No afecta nada!!
    
    A_new_base(2) = Constante - Rx / (distance * speed * Ry );
%    A_new_base(2) = - Rx / (distance * speed * Ry );

end

A_new_base(3) = 0;

% Ahora volvemos a la base inicial.

A(1) = A_new_base(1) * v_norm(1) + A_new_base(2) * e_2(1);

A(2) = A_new_base(1) * v_norm(2) + A_new_base(2) * e_2(2);

A(3) = A_new_base(1) * v_norm(3) + A_new_base(2) * e_2(3);
