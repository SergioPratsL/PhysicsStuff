
function Y = fSigma(V)
% Usar velocida normalizada!

veloc = norm(V);

Y = 1 / sqrt( 1 - veloc^2);


