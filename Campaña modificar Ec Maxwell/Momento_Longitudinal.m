function [ p_long ] = Momento_Longitudinal( TMaxwell )
% Funcion que diagonaliza el tensor de Maxwell y obtiene el momento 
%a partir del autovector  y autovalor de la componente positiva 

[V, D] = eig(TMaxwell);

if D(1,1) > 0
    p_long = V(1:3, 1) * D(1,1);
elseif D(2,2) > 0
    p_long = V(1:3, 2) * D(2,2);
elseif D(3,3) > 0
    p_long = V(1:3, 3) * D(3,3);    
else
    p_long = [0; 0; 0]
end

p_long = p_long';

