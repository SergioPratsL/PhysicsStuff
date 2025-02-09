function [ V, A ] = Potencial_LW( R, v )
% R es la distancia retardada y v la velocidad de la carga que genera el
% potencial, con velocidad de la luz normalizada a 1 y pasando de la carga
% y otras constantes

V = 1 /(norm(R) - dot(R,v));

A = v  / (norm(R) - dot(R,v));

end

