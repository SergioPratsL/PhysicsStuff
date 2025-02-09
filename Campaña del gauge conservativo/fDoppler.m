function Dop = fDoppler(v, R)
% v es la velocidad con la que ves que se mueve el emisor
% R es la distancia retardada de la fuente en el momento de calcular el
% doppler
% v y R son ambos vectores de 3 dimensiones [x, y, z]
% Velocidad como siempre normalizada c=1

veloc = norm(v);

Sigma = fGamma(veloc);

R_norm = R / norm(R);

Rv = dot(R_norm, v);

Dop = 1/(Sigma) * 1/(1+Rv);



