function [pendiente_new, dif_p_x, cambiaPsi] = MomentoRebote( pendiente, tipo, tan_alfa, Psi )
%[Rho_new, dif_p_x] = MomentoRebote( alfa, pendiente )


%if Psi <= pi/2
%    A = [1, pendiente];
%else
%    A = [-1, -pendiente];
%end
% El efecto de  Psi se ve reflejado en "alfa_efectiva"
A = [1, pendiente];

A_norm = A / norm(A);
    

if tipo == 'rampa'

    alfa_efectiva = atan( tan_alfa * cos(Psi) );
    
% Rotar al eje de alfa (signo rotacion estaba mal)
    A_rot(1) = cos(alfa_efectiva) * A_norm(1) + sin(alfa_efectiva) * A_norm(2);
    A_rot(2) = cos(alfa_efectiva) * A_norm(2) - sin(alfa_efectiva) * A_norm(1);

% Invertir direccion y'
    A_rot(2) = - A_rot(2);

% Volver al eje original
    A_new(1) = cos(alfa_efectiva) * A_rot(1) - sin(alfa_efectiva) * A_rot(2);
    A_new(2) = cos(alfa_efectiva) * A_rot(2) + sin(alfa_efectiva) * A_rot(1);

% Momento transferido al suelo
    dif_p_x = A_norm(1) - A_new(1);
    pendiente_new = A_new(2) / A_new(1);
    
    if sign(A_new(1)) ~= sign(A_norm(1))
        cambiaPsi = 1;
    else
        cambiaPsi = 0;
    end
    
elseif tipo == 'muro_'     
    dif_p_x = - 2 * A_norm(1);
    pendiente_new = - pendiente;   
    
end