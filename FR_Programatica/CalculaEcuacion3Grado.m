function x = CalculaEcuacion3Grado(b,c)
% Resuelve x + b*x^3 = c

% 01.08.2022.

% fact1 = ( (3*b)^(1/2)*((27*c^2+4)^(1/2) + 9*b*c))^(1/3);
% term_1 =  fact1 / (2^(1/3) * (3*b)^(2/3));
% term_2 = (2/(3*b))^(1/3) / fact1;
% x = term_1 - term_2;

%Copiado de internet, va mal.

% Bolzano, siempre serán valores positivos.

max_err = c / 20;
max_iter = 100;

dif_percent = 0.25;

n = 0;
val_ini = c;

%%% CANCELADO, LA FUERZA DE REACCIÓN HA CAIDO, ¡¡¡¡HA CAIDOOOO!!!!
%%% AL FINAL SÓLO HIZO FALTA CAMBIAR EL PUTO SÍGNO DEL a^2*a XDDDDD

% Para muy altas energías quizá haría falta pero de momento lo dejo aquí

end

