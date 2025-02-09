function dv_red = DifVelRed( dv, v_vect, cad );
% Package: Prueba de la aceleracion dual

% Obtiene la contribucion que un acelerin hara en la velocidad total
% aplicandole los sigmas pertinentes en la base de dicha velocidad

T = CambioBaseEnh(dv, v_vect);

T_inv = T^(-1);


v = norm(v_vect);
Sigma = 1 / sqrt(1-v^2);

if cad == 'E'           % Emisor
    M = - [1/Sigma^2,0, 0; 0, 1, 0; 0, 0, 1];       % Tomasito...
elseif cad == 'R'       % Receptor
    M = [1/Sigma^2,0, 0; 0, 1/Sigma, 0; 0, 0, 1/Sigma];
end 
    
% dv_new_base_red = (dv * T) * M;
% 
% dv_red = T_inv * dv_new_base_red';
% 
% dv_red = dv_red';

% Intento 2
% dv_red = T_inv * (M * (dv * T)');


% Desesperao...
p = (T * dv' )';
p2 = p * M;
dv_red = p2 * T_inv';



aux = 13;


