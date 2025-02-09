
% Prueba 1. Finalmente OK
Rx = 0;
Ry = 1;
v = 0.4;
a = [1, 0, 0];

% Prueba 2. OK
Rx = 0;
Ry = 1;
v = 0.4;
a = [0, 0, 1];

% Prueba 3. Por fin!!
Rx = 0;
Ry = 1;
v = 0.4;
a = [0, 1, 0];


% Prueba 4. OK
Rx = -1.1;
Ry = 1;
v = 0.4;
a = [1, 0, 0];


% Prueba 5. OK
Rx = -1.1;
Ry = 1.6;
v = -0.4;
a = [0, 1, 0];


% Prueba 6. OK
Rx = -2.1;
Ry = -1.6;
v = -0.21;
a = [0, 0, 1];


% Prueba 7. (Cambio signo de la aceleracion). OK
Rx = -1.1;
Ry = 1;
v = 0.4;
a = [-1, 0, 0];


% Prueba 8. (Cambio signo a). OK
Rx = -1.1;
Ry = 1.6;
v = 0.4;
a = [0, -1, 0];


% Prueba 9. (Cambio signo a). OK
Rx = -2.1;
Ry = -1.6;
v = -0.26;
a = [0, 0, -1];

% PRUEBAS "GENERALES". 

% Prueba 10. OK
Rx = 1.4;
Ry = -0.9;
v = -0.43;
a = [2.1, 0, -1];

% Prueba 11. OK
Rx = 2.4;
Ry = -3.9;
v = -0.33;
a = [0.4, -0.6, 0.8];

% Prueba 12. OK
Rx = 0;
Ry = -2;
v = -0.33;
a = [0.5, -1.2, 1];

% Prueba 13. OK
Rx = 0;
Ry = -1;
v = -0.33;
a = [0, 1, 1];

% Prueba 14. OK
Rx = 4;
Ry = -2;
v = -0.33;
a = [0.5, 0, 1];


% Pruebas desde la lejania (X >> Y)

% Prueba 15. OK, el error crece con la lejania
Rx = 1;
Ry = 0.0001;
v = 0.37;
a = [1, 0, 0];

% Prueba 16. % El error crece al cuadrado que X/Y!! Muy preocupante
Rx = 1;
Ry = 0.0001;
v = 0.37;
a = [0, -1, 0];


% Prueba 17. % El error crece al cuadrado que X/Y!! Muy preocupante
Rx = 1;
Ry = 0.0001;
v = 0.37;
a = [0, 0, 1];


% Baja velocidad
% Prueba 18. El error crece con el descenso de la velocidad al cuadrado!
Rx = 1;
Ry = 1;
v = 0.0001;
a = [1, 0, 0];


% Prueba 19. El error crece con la inversa a la velocidad (pongo ceros)
Rx = 1;
Ry = 1;
v = 0.0001;
a = [0, 1, 0];


% Prueba 20. El error crece con el descenso de la velocidad al cuadrado!
Rx = 1;
Ry = 1;
v = 0.00001;
a = [0, 0, -1];


%Alta velocidad
% Prueba 21. % El error decrece aun mas con mas 9s
Rx = 1;
Ry = 1;
v = 0.9999;
a = [1, 0, 0];


% Prueba 22. % El error crece proporcionalmente a los 9s que hay
Rx = 1;
Ry = 1;
v = 0.9999;
a = [0, 1, 0];


% Prueba 23. % El error crece proporcionalmente a los 9s que hay
Rx = 1;
Ry = 1;
v = 0.999;
a = [0, 0, 1];


v_0 = [0, 0, 0];

R = [Rx, Ry, 0];

R_fut2 = PostRetFut( R, v, 'f');

[E_deriv_LW, B_deriv_LW] = CampoRadiadoYinn( Rx, Ry, v, a)

% 16.12.2015. Usa la distancia futura so SUBNORMAL!!
%[E, B] = CampoRadiadoCargaAcelerada(R, v_0, a)
[E, B] = CampoRadiadoCargaAcelerada(R_fut2, v_0, a)

Dif_E = E - E_deriv_LW 

Dif_B = B - B_deriv_LW