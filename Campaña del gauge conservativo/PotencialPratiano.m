function [V, A] = PotencialPratiano(R, v)

% 22.11.2015. Este potencial es malo, para el bueno,
% ver PotencialDeYinn :P.

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
    
    if v_norm(2) ~= 0  ||  v_norm(3) ~= 0
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

if Rx == distance * speed
    A_new_base(2) = 0;
    
else   
    A_new_base(2) = - Ry / distance^3 * (distance * Rx * (1-speed^2) -speed * Ry^2) / (Rx - distance * speed) ;

end

A_new_base(3) = 0,

% Ahora volvemos a la base inicial.

A(1) = A_new_base(1) * v_norm(1) + A_new_base(2) * e_2(1);

A(2) = A_new_base(1) * v_norm(2) + A_new_base(2) * e_2(2);

A(3) = A_new_base(1) * v_norm(3) + A_new_base(2) * e_2(3);
