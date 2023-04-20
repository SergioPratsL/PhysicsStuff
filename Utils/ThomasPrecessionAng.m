function [dir_rot, angulo] = ThomasPrecessionAng( Sigma, dv, v )
% Sigma está asociado a v, se puede pasar para ahorrar operaciones.

if Sigma == 0
    Sigma = fSigma(v);
end

v_rot = Sigma^2 /(Sigma + 1) * cross(dv, v);

angulo = norm(v_rot);

if angulo == 0
    dir_rot = [1,0,0];  % Irrelevante
else
    dir_rot = v_rot / angulo;    
end


