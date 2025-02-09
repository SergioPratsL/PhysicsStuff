function TMomento_Energia = Tens_Momento_Energia(E, B)
% Función que calcula el tensor de Mawell.
% Se omiten todas las constantes para trabajar de forma simplificada

M = (E' * E) + (B' * B);

if ndims(M) == 1
    printf 'Algo ha ido mal, los vectores no se pasaron correctamente'
end

cof = (1/2) * (norm(E)^2 + norm(B)^2);

M(1,1) = M(1,1) - cof;
M(2,2) = M(2,2) - cof;
M(3,3) = M(3,3) - cof;

u = norm(E)^2 * 0.5 + norm(B)^2 * 0.5;

S = cross(E,B);

T = zeros(4,4);
T(1,1) = u;
T(2:4, 1) = S;
T(1, 2:4) = S';
T(2:4, 2:4) = - M;

TMomento_Energia = T;

