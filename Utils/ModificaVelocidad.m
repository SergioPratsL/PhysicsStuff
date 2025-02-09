function [v_new, Sigma_new] = ModificaVelocidad( v, dp );
% v es la velocidad inicial y dp es a*dt, la variación de momento vista
% desde el laboratorio, claro.
% obviamente hago todo con masa en reposo = 1.

Sigma = fGamma(v);

% Esto va bien en línea recta si no, no.
% p = Sigma * v;
% 
% p_new = p + dp;
% 
% Sigma_new = sqrt(1 + p_new(1)^2 + p_new(2)^2 + p_new(3)^2);
% 
% v_new = p_new / Sigma_new;

% No mejoró nada sino que dio el mismo resultado.
% El boost es a_prop * dt / Sigma, lo mismo que dp/Sigma
% dv = dp/Sigma;
% % El menos porque el lab se mueve a -v y esto es un boost en el que acelera
% v_ini = [1, -v];
% 
% v_u_new = Boost(v_ini, dv);
% v_new = -v_u_new(2:4) / v_u_new(1);
% 
% Sigma_new = fGamma(v_new);

% Cuánto daño hace la precesión de Thomas.
% if (norm(cross(v,v_new)) > 10^-7)
%     rotMatrix = GetThomasRotMatrix(v, v_new);
%     v_new = v_new * rotMatrix';
%     %v_new = PutoThomas(v, v_new);
% end
    

% Una tercera forma de hacerlo que no debería requerir rotación
% Boost de A_new a A, luego al lab.

dv = dp / Sigma;

dv_u = [1, dv]; % Sigma debería ser despreciable
v_u = Boost(dv_u, - v);

v_new = v_u(2:4) / v_u(1);

Sigma_new = fGamma(v_new);