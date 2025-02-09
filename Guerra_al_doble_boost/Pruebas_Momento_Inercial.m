clear;

vA = [0.8, 0, 0];
%vC = [-0.001, 0, 0];
vC = [0, -0.001, 0];

SigmaA = fGamma(vA);
pA = vA * SigmaA

vAC = Vel_Addition_Law( vA, vC )

SigmaAC = fGamma(vAC)
pAC = vAC * SigmaAC

% Sigma veces mayor tanto con boost paralelo como perpendicular
dif_p = pAC - pA       

