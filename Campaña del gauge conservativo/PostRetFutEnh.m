function R_out = PostRetFutEnh(R_vect, v_vect, cad)
% Package: prueba geodesica
% Este script te da la posicion futura o pasada a partir de la posicion
% actual para una velocidad v que no tiene que ser en direccion X ni R
% tiene porque estar en el plano XY

% La idea es simple, cambiar la base para que v_vect sea X.

v = norm(v_vect);
R = norm(R_vect);

% Excepciones que no hay que tratar.
if v == 0
    R_out = R_vect;
    return
elseif R == 0
    R_out = R_vect
    return
end

e1 = v_vect / v;

% Eleccion de e2, debe estar en el plano R_vect, v_vect

R_vect_aux = R_vect - dot(R_vect, e1) * e1;
R_aux = norm(R_vect_aux);

if norm(R_vect_aux) ~= 0
    
    e2 = R_vect_aux / R_aux;
    
% R y v alineadas, elegir cualquier plano    
else
    R_vect_aux2 = [0,1,0];
    if dot(e1, R_aux) ~= 1
        R_vect_aux = R_vect_aux2 - dot(e1, R_aux) * e1;
    else        
        R_vect_aux2 = [1,0,0];
        R_vect_aux = R_vect_aux2 - dot(e1, R_aux) * e1;
    end
        
    e2 = R_vect_aux / norm(R_vect_aux);
end

e3 = cross(e1, e2);
e3 = e3 / norm(e3);        % Por si las moscas

% Expresamos los vectores en la base privilegiada
v_new_base = [v, 0, 0];

R1 = dot(R_vect, e1);
R2 = dot(R_vect, e2);
R3 = dot(R_vect, e3);      % Deberia ser 0!


R_new_base = [R1, R2, R3];


R_out_new_base = PostRetFut(R_new_base, v, cad);

R_out_1 = R_out_new_base(1) * e1(1) + R_out_new_base(2) * e2(1) + R_out_new_base(3) * e3(1);
R_out_2 = R_out_new_base(1) * e1(2) + R_out_new_base(2) * e2(2) + R_out_new_base(3) * e3(2);
R_out_3 = R_out_new_base(1) * e1(3) + R_out_new_base(2) * e2(3) + R_out_new_base(3) * e3(3);


R_out = [R_out_1, R_out_2, R_out_3];


