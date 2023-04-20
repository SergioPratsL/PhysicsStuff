

function [Eo, Bo] = Boost_EM(Ei, Bi, V)
% Transforma el campo EM entre dos SRIs
% Velocidad normalizada! (c=1)
% Unidades de E y B omitidas!

% Hay que probar que la rotacion caiga sobre el eje Z ya que 
% los ejes XY se unen contra Z y asi se puede perder la informacion.
% Se tendran que repetir las pruebas

% Eo = [Ex, Ey, Ez]
% Bo = [Bx, By, Bz]
% Ei = [Ex, Ey, Ez]
% Bi = [Bx, By, Bz]
% V = [vx, vy, vz]

dims = size(Ei);

if dims(1) ~= 1  || dims(2) ~= 3
    disp('mal')
    return
end     

dims = size(Bi);

if dims(1) ~= 1  || dims(2) ~= 3
    disp('mal')
    return
end     

dims = size(V);

if dims(1) ~= 1 || dims(2) ~= 3 
    disp('mal')
    return
end 
    
if norm(V) >= 1
    disp('mal')
    return
end 

veloc = norm(V);
%Sigma = 1 / sqrt( 1 - veloc^2);
Sigma = fSigma(veloc);

if veloc == 0
    Eo = Ei;
    Bo = Bi;
    return 
end 

% Saco las componentes del vector
e3 = V / veloc;
[e1, e2] = gPerpendicular( e3);



% Obtenemos el producto vectorial, ya multiplicado por Sigma
ExV = Sigma * cross( V, Ei );
BxV = Sigma * cross( V, Bi );

% Aplicar la Sigma a las componentes perpendiculares
E_new_base = CambioBase(e1, e2, e3, Ei);
E_new_base(1) = Sigma * E_new_base(1);
E_new_base(2) = Sigma * E_new_base(2);
% Volver a la base normal
E_Sigma = E_new_base(1)*e1 + E_new_base(2)*e2  + E_new_base(3)*e3;

% Ahora lo mismo para el campo magnetico
B_new_base = CambioBase(e1, e2, e3, Bi);
B_new_base(1) = Sigma * B_new_base(1);
B_new_base(2) = Sigma * B_new_base(2);
% Volver a la base normal
B_Sigma = B_new_base(1)*e1 + B_new_base(2)*e2  + B_new_base(3)*e3;


% Finalmente, obtener los vectores buenos
Eo = E_Sigma + BxV;
Bo = B_Sigma - ExV;


% Blasfemias varias para lograr que vaya la comparativa del campo radiado
%Eo = Eo * norm(E_Sigma) / norm(Eo);
%Bo = Bo * norm(B_Sigma) / norm(Bo);
