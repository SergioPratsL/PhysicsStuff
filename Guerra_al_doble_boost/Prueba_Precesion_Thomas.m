clear;

% Quiero verificar que la formula de la precesión de Thomas va bien

% Velocidad inicial de la particula respecto al laboratorio (L)
%vA = [0.4, 0, 0];
% Variación de velocidad en el SRI local de la partícula
%dv = [0, 0.001, 0];

% vA = [0.6, 0, 0];
% dv = [-0.002, 0.001, 0];

% vA = [0.6, 0, 0];
% dv = [-0.002, 0.001, 0.0014];

vA = [0.95, 0, 0];
dv = [-0.002, 0.001, 0.0014];


vLabA_new = Vel_Addition_Law( dv, -vA )
vALab_new = Vel_Addition_Law( -vA, dv )

dv_lab = vLabA_new - vA;

vALab_new_rot = ThomasPrecession(0, dv_lab, vA, vALab_new)
