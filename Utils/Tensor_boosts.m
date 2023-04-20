function T_boost_v = Tensor_boosts(v)
% Da el tensor con el que puedes transformar boosts, se le pasa solo una
% velocidad (asume c=1)
% Hace boosts de tu propio sistema a otro sistema que se mueve a velocidad
% v, para hacer el boost en sentido contrario trasponer la matriz!

% Fuente: https://en.wikipedia.org/wiki/Lorentz_transformation#Contravariant_vectors

B = zeros(4);

speed = norm(v);
v_norm = v / speed;
Sigma =  1 / sqrt(1-speed^2);

B(1,1) = Sigma;

B(2,1) =  - Sigma * v(1);
B(3,1) =  - Sigma * v(2);
B(4,1) =  - Sigma * v(3);

B(1,2) =  - Sigma * v(1);
B(1,3) =  - Sigma * v(2);
B(1,4) =  - Sigma * v(3);

B(2,2) = 1 + (Sigma - 1) * v_norm(1)^2;
B(3,2) = (Sigma - 1) * v_norm(1)* v_norm(2);
B(4,2) = (Sigma - 1) * v_norm(1)* v_norm(3);

B(2,3) = (Sigma - 1) * v_norm(2)* v_norm(1);
B(3,3) = 1 + (Sigma - 1) * v_norm(2)^2;
B(4,3) = (Sigma - 1) * v_norm(2)* v_norm(3);


B(2,4) = (Sigma - 1) * v_norm(3)* v_norm(1);
B(3,4) = (Sigma - 1) * v_norm(3)* v_norm(2);
B(4,4) = 1 + (Sigma - 1) * v_norm(3)^2;

T_boost_v = B;
