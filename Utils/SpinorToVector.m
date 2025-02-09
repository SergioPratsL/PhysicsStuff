function V = SpinorToVector(bs)
% https://en.wikipedia.org/wiki/Pauli_matrices#Pauli_vectors

% A sería una matriz suma de matrices de Pauli por un vector
% A = [1, 0; 0, -1];
% V daría los autovalores (-1, 1) mientras que D serían los bispinors...
% [V, D] = eig(A)

% Muy bien, pero quiero lo contrario...
% Una matriz tal que V*M = 1

bs = bs / norm(bs);

sx = MatrizPauli(1);
sy = MatrizPauli(2);
sz = MatrizPauli(3);

x = bs * sx * bs';
y = bs * sy * bs';
z = bs * sz * bs';

V = [x, y, z];


% Pruebas
%A = SpinorToVector([i,0])    [0, 0, 1]
%A = SpinorToVector([i,3])    [0    0.6000   -0.8000]
% A = SpinorToVector([i,1.08i]) [0.9970         0   -0.0768]
% A = SpinorToVector([(0.7070 + 0.01*i), (-0.01 - 0.707*i) ])     [-0.0283    0.9996         0]
