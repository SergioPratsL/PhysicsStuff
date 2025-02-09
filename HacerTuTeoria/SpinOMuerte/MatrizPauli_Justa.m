% Astracanada
function Sigma = MatrizPauli_Justa(index)

% Matriz 4x4, causará un error, ocurre si index es incorrecto
Sigma = [1,1,1,1;1,1,1,1;1,1,1,1]; 

if index == 1
    Sigma = [0,1+1i; 1+1i,0]/sqrt(2);
end

if index == 2
    Sigma = [0, -1+1i; 1-1i, 0]/sqrt(2);
end

if index == 3
    Sigma = [1,0; 0, -1];
end
