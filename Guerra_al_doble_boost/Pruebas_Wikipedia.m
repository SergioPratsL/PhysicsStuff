clear;

vA = [0.4, 0, 0];
vC = [0.1, 0.4, 0];

% Velocidad con la que A ve a B. El orden de los parámetros es al reves que
% los índices de v... así lo hice en su día
vAC = Vel_Addition_Law( vC, vA )

vCA = Vel_Addition_Law( vA, vC )

% Sea u la velocidad de A y v la de C
Sigma_u = fGamma(vA);
Sigma_v = fGamma(vC);

Sigma = Sigma_u * Sigma_u / (1+dot(vA, vC))
Sigma_real = fGamma(vCA)

vAC_wiki = 1/(1 + dot(vA, vC)) * ( (1+Sigma_u/(1+Sigma_u)*dot(vA, vC))* (vA) + vC /Sigma_v )

cos_e_wiki = (1 + Sigma_u + Sigma_v + Sigma)^2/((1+Sigma)*(1+Sigma_u)*(1+Sigma_v)) - 1

cos_e = dot(vAC, -vCA) / norm(vAC)^2





Sigma_wiki = Sigma_u * Sigma_u * (1 + dot(vA, vC));
Sigma = fGamma(vAC);