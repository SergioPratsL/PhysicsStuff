function [Ei, Bi, Er, Br] = Campo_EM_Total_No_Units(R, v, a, t_unit)
% 2021 y parece mentira que en toda la serie de Matlabs que se remonta a
% 2014 no haya una puta función que saque el campo radiado :P.

% t_unit es el tiempo unitario para la velocidad, se entiende que 
% la unidad de distancia es el tiempo que hace la luz en t_unit y la
% velocidad de la luz es 1 (1 * d_unit/t_unit).
% Aclarado esto, en verdad t_unit no se usa para nada aquí.

[Ei, Bi] = CampoInducido_sin_unidades(R, v);

[Er, Br] = CampoRadiadoCargaAcelerada_No_Units(R, v, a);

end


