function A = PotencialPratiano(R, v)

% Funcion que devuelve el potencial pratiano en el SRI propio de la carga
% sobre la que evaluas en potencial

% R es el vector que va desde el origen al destino (el destino es el punto
% que evaluamos.


% v es un numero (puede ser negativo), puesto que en esta version se asume
% que la velocidad es sobre el eje X (el trabajo de cambiar de base no se
% hace dentro de esta funcion.

% La velocidad debe estar normalizada respecto la velocidad de la luz!

% Lo primero sacar el campo inducido, para lo que hace falta la posicion
% retardada!

Rr = DRetardada(R, v);

[E, B] = CampoInducido_sin_unidades(Rr, v);

% Doppler
Doppler = 1 / (1 - Rr(1)*v/norm(Rr) );
dv = [0, Doppler, 0];

[E_rad, B_rad] = CampoRadiadoCargaAcelerada( Rr, v, dv)

% El campo electrico depende de la distancia retardada pero las R's
% que salen en el potencial dependen de la actual!!

A(1) = - R(1) / v * E(1) - R(2) / v * E(2) + R(2) * B(3);


A_21 = ( R(2) - R(1) ) / v * E(1);

% Pendiente confirmar que con Rx = 0 no petara esta parte...
A_22 = R(2) / R(1) * (v * Erad(2) + R(2) * B(3));

A(2) = A_21 + A22;

A(3) = 0;

% Petara por todas partes...