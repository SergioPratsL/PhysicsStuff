function Gamma = MatrizGamma(index)

% Matriz 3x3, causará un error, ocurre si index es incorrecto
Gamma = [1,1,1;1,1,1;1,1,1]; 

if index == 0
    Gamma = [1,0,0,0; 0,1,0,0; 0,0,-1,0; 0,0,0,-1];
end

if index == 1
    Gamma = [0,0,0,1; 0,0,1,0; 0,-1,0,0; -1,0,0,0];
end

if index == 2
    Gamma = [0,0,0,-1i; 0,0,1i,0; 0,1i,0,0; -1i,0,0,0];
end

if index == 3
    Gamma = [0,0,1,0; 0,0,0,-1; -1,0,0,0; 0,1,0,0];
end


if index == 5
    Gamma = [0,0,1,0; 0,0,0,1; 1,0,0,0; 0,1,0,0];
end
