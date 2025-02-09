function Sigma = MatrizPauli(index)
% https://en.wikipedia.org/wiki/Pauli_matrices

% Matriz 4x4, causará un error, ocurre si index es incorrecto
Sigma = [1,1,1,1;1,1,1,1;1,1,1,1]; 

if index == 1
    Sigma = [0,1; 1,0];
end

if index == 2
    Sigma = [0, -1i; 1i, 0];
end

if index == 3
    Sigma = [1,0; 0, -1];
end
