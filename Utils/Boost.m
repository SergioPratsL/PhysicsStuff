
function R2 = Boost(R1, V)
% R1 = [T, X, Y, Z]
% V = [vx, vy, vz]
% R2 = [T, X, Y, Z]

% V es la velocidad con la que veo que se mueve el cuerpo para el cual
% quiero transformar el vector R1, que esta sacado de mi sistema

dims = size(R1);

if dims(1) ~= 1  || dims(2) ~= 4
    disp('mal') 
    return
end     
    
dims = size(V);

if dims(1) ~= 1 || dims(2) ~= 3 
    disp('mal')
    return
end 

veloc = norm(V);

if veloc >= 1
    disp('mal')
    return
end 
    
if veloc == 0
    R2 = R1;
    return 
end 


Sigma = 1 / sqrt( 1 - veloc^2);

% Obtengo la direccion de la velocidad y sus ortogonales
e3 = V / veloc;
[e1, e2] = gPerpendicular( e3);

R1_esp = [R1(2)  R1(3)  R1(4)];

R1_new_base = CambioBase(e1, e2, e3, R1_esp);

% Inicializo R2 y la parte espacial
R2 = [0, 0, 0, 0];
R2_new_base = [0, 0, 0];

% Tiempo en el otro SRI
R2(1) = Sigma * (R1(1) - veloc * R1_new_base(3));

% Componentes espaciales de la base privilegiada
R2_new_base(1) = R1_new_base(1);
R2_new_base(2) = R1_new_base(2);
R2_new_base(3) = Sigma * (R1_new_base(3) - veloc * R1(1));

R2_esp = R2_new_base(1) * e1 + R2_new_base(2) * e2 + R2_new_base(3) * e3;

R2(2) = R2_esp(1);
R2(3) = R2_esp(2);
R2(4) = R2_esp(3);
