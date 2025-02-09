
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
Sigma = 1 / sqrt( 1 - veloc^2);

if veloc == 0
    Eo = Ei
    Bo = Bi
    return 
end 

% Sacamos el angulo de la velocidad respecto a x.
veloc = norm(V);
V_norm = V / veloc;

if V(2) == 0 && V(1) >= 0
    Fi = 0;
elseif V(2) == 0 && V(1) < 0
    Fi = pi / 2;
else
% Mismo angulo si x < 0 y y < 0 que si x > 0 y y > 0, me parece bien 
% pues es necesaria la misma rotacion para alinearse con un eje
    Fi = atan( V(2) / V(1) );
    
end 
    
XY = [V(1)  V(2)];

norm_XY = norm(XY);

if norm_XY == 0 
    if V(3) >= 0
        Theta = pi/2;
    else
        Theta = -pi/2;
    end 
else     
    
    Theta = atan( V(3) / norm_XY );

end

if Theta > 0
    dTheta = pi/2 - Theta
else
    dTheta = -pi/2 - Theta
end

% En polares no hay que usar la formula de Rodrigues! 
% hay que expresar los campos en polares y rotarlos.

% Rotamos la velocidad (debe quedar en direccion x pero quiero comprobarlo
% V_rot = rotacion_polar( V, -Fi, -Theta );   (old)
V_rot = rotacion_polar( V, 0, dTheta);

% Rotamos los campos para que el boost sea en direccion Z.
% Ei_rot = rotacion_polar( Ei, -Fi, -Theta ); (old)
% Bi_rot = rotacion_polar( Bi, -Fi, -Theta ); (old)
Ei_rot = rotacion_polar( Ei, 0, dTheta );
Bi_rot = rotacion_polar( Bi, 0, dTheta );

% Obtenemos el producto vectorial
ExV = cross( V_rot, Ei_rot );
BxV = cross( V_rot, Bi_rot );

% Obtenemos los vectores transformados (y rotados)
% Eo_rot = [Ei_rot(1),  Sigma * ( Ei_rot(2) + BxV(2) ), Sigma * ( Ei_rot(3) + BxV(3) ) ];
% Bo_rot = [Bi_rot(1),  Sigma * ( Bi_rot(2) - ExV(2) ), Sigma * ( Bi_rot(3) - ExV(3) ) ]; 

Eo_rot = [Sigma * ( Ei_rot(1) + BxV(1) ),  Sigma * ( Ei_rot(2) + BxV(2) ), Ei_rot(3) ];
Bo_rot = [Sigma * ( Bi_rot(1) - ExV(1) ),  Sigma * ( Bi_rot(2) - ExV(2) ), Bi_rot(3) ];

% Volvemos a girar los vectores a su eje original
%Eo = rotacion_polar( Eo_rot, Fi, Theta); (old)
%Bo = rotacion_polar( Bo_rot, Fi, Theta); (old)

Eo = rotacion_polar( Eo_rot, 0, -dTheta);
Bo = rotacion_polar( Bo_rot, 0, -dTheta);

