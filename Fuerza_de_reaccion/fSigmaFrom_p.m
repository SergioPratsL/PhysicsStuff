function Sigma = fGammaFrom_p(p)
% El objetivo de esto es sacar Sigma pero no de la velocidad sino del 
% momento, asumiendo que p está calculado sobre una masa en reposo = 1.
% c = 1, como de costumbre

p_norm = norm(p);

v_norm = sqrt(p_norm^2/(1+p_norm^2));

Sigma = 1 / sqrt(1 - v_norm^2);