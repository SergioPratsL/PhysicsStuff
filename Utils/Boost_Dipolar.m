

function [pf, mf] = Boost_Dipolar(pi, mi, v)
% Adaptado de Boost_EM

dims = size(pi);

if dims(1) ~= 1  || dims(2) ~= 3
    disp('mal')
    return
end     

dims = size(mi);

if dims(1) ~= 1  || dims(2) ~= 3
    disp('mal')
    return
end     

dims = size(v);

if dims(1) ~= 1 || dims(2) ~= 3 
    disp('mal')
    return
end 
    
if norm(v) >= 1
    disp('mal')
    return
end 

veloc = norm(v);
%Sigma = 1 / sqrt( 1 - veloc^2);
Sigma = fGamma(veloc);

if veloc == 0
    pf = pi;
    mf = mi;
    return 
end 

% Saco las componentes del vector
e3 = v / veloc;
[e1, e2] = gPerpendicular(e3);


% Obtenemos el producto vectorial, ya multiplicado por Sigma
vxp = Sigma * cross( v, pi );
vxm = Sigma * cross( v, mi );

% Aplicar la Sigma a las componentes perpendiculares (e3 es la dir de la
% velocidad)
p_new_base = CambioBase(e1, e2, e3, pi);
p_new_base(1) = Sigma * p_new_base(1);
p_new_base(2) = Sigma * p_new_base(2);
% Volver a la base normal
p_Sigma = p_new_base(1)*e1 + p_new_base(2)*e2  + p_new_base(3)*e3;

% Ahora lo mismo para el campo magnetico
m_new_base = CambioBase(e1, e2, e3, mi);
m_new_base(1) = Sigma * m_new_base(1);
m_new_base(2) = Sigma * m_new_base(2);
% Volver a la base normal
m_Sigma = m_new_base(1)*e1 + m_new_base(2)*e2  + m_new_base(3)*e3;


% Finalmente, obtener los vectores buenos

% Algunos signos invertidos.
pf = p_Sigma + vxm;
mf = m_Sigma - vxp;

% Muy malo :(
% pf = p_Sigma - vxm;
% mf = m_Sigma + vxp;

% Muy malo :(
% pf = p_Sigma - vxm;
% mf = m_Sigma - vxp;


% Signos cambiado máss otro que está muy mal.
%pf = p_Sigma + vxm;
%mf = m_Sigma + vxp;
