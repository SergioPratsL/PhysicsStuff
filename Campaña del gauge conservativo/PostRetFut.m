function R_out = PostRetFut(R, v, cad)

% Este script te da la posicion futura o pasada a partir de la posicion
% actual.
% Orientado a las pruebas finales de la teoria del potencial Unico (Gauge
% Killer).

% Premisas: v es un numero con la velocidad en direccion x normalizada a la
% velocidad de la luz.

% R es un vector de 3 valores pero Z debe valer 0!
% R cuando se saca la velocidad pasada va del cuerpo que se mueve al que
% esta quieto mientras que para sacar la velocidad futura va del cuerpo
% quieto hacia el que se mueve.

R_out = [0, 0, 0];

R_out(2) = R(2);
R_out(3) = 0;

Sigma = 1 / sqrt(1 - v^2);

Raiz = Sigma * v * sqrt((Sigma^2 * R(1)^2 + R(2)^2 ));

Parte_1 = Sigma^2 * R(1);

vel = [v, 0, 0];

if cad == 'r'
    R_out(1) = Parte_1 - Raiz;

    Rr_verif = R - norm(R_out) * vel;
    dif_verif = R_out - Rr_verif;
    
    
elseif cad == 'f'
    R_out(1) = Parte_1 + Raiz;
    
    Rr_verif = R + norm(R_out) * vel;
    dif_verif = R_out - Rr_verif;
    
end


