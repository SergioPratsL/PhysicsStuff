
function TMaxwell = TensorDeMaxwell(E, B)
% Función que calcula el tensor de Mawell.
% Se omiten todas las constantes para trabajar de forma simplificada

%M = E' * B;        % Argg!

M = (E' * E) + (B' * B);

if ndims(M) == 1
    printf 'Algo ha ido mal, los vectores no se pasaron correctamente'
end

cof = (1/2) * (norm(E)^2 + norm(B)^2);

M(1,1) = M(1,1) - cof;
M(2,2) = M(2,2) - cof;
M(3,3) = M(3,3) - cof;

TMaxwell = M;


