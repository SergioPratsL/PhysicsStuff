function T = CambioBaseEnh(R_vect, v_vect)
% Package: aceleracion dual

% Define una transformacion para que v_vect quede en direccion X
% v_vect pasa a ser +X

v = norm(v_vect);
v_norm = v_vect / v;

if norm(v) == 0
    return
end 

R = norm(R_vect);
R_norm = R_vect / R;

if norm(R) == 0        % :@
    return
end

e_1 = v_norm;


aux = cross(R_vect, v_vect);

if norm(aux) ~= 0
% if aux ~= 0
    e_2_not_norm = R_vect - v_norm * dot(R_vect, v_norm);
else
    v_aux = [0, 1, 0];
    aux2 = cross( v_vect, v_aux);
    
   if norm(aux2) ~= 0
%    if aux2 ~= 0
        e_2_not_norm = cross(v_norm, v_aux);
    else
        v_aux2 = [0, 0, 1];
        e_2_not_norm = cross(v_norm, v_aux2);
    end
end
    
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
    
    v_norm_ort = vector_ortogonal / norm_v_ort;
    
    vect_ortogonal = cross( v_norm_ort, vect_test );
    
    vector_ortonormal = vector_ortogonal / norm(vector_ortogonal);    
    
end

if norm(e_2_not_norm) == 0
    e_2 = cross( v_norm, vector_ortonormal );
    e_2 = e_2 / norm(e_2);      % Por si las moscas
else
    e_2 = e_2_not_norm / norm(e_2_not_norm);
end

e_3 = cross(e_1, e_2);


T = [e_1; e_2; e_3];

