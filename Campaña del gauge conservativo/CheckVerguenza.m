function result = CheckVerguenza( Rx, Ry, v )

% Uso de la fuerza bruta para una expresion que no se transformar...
% La velocidad siempre en X, por lo que es un numero, y esta normalizado
% respecto a la velocidad de la luz.

R = sqrt( Rx^2 + Ry^2 );

Sigma2 = 1 / (1 - v^2);

Fact1 = Ry / (Sigma2 * R^5 * (Rx - R*v)^2);


Sum1 = ( 3 * R * Rx^2 - 3 * Sigma2 * v * Rx * Ry^2 ) * (Rx - R*v);

Sum2 = - R * Rx^3 + v * R^4 - Sigma2 * v * R^2 * Ry^2 + Sigma2 * v^2 * R * Rx * Ry^2;


result_ideal = - Ry / (R^3*v)

result = Fact1 * (Sum1 + Sum2);


