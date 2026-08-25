
gt = MatrizGamma(0);
g5 = MatrizGamma(5);
[at, ax, ay, az] = MatricesAlfa();

Sz = MatrizSpin_4_4(3);


%Pruebas:
%v = [0, 0, 0];
spinor_base = [1, 0]; 
%spinor_base = [1, 1]; 
%spinor_base = [cos(pi/8), sin(pi/8)]    


%v = [0, 0, 0.4];
%v = [0.4, 0, 0];
%v = -[0.3, 0, -0.350];
v = -[0.3, 0.28, -0.350];

p = fGamma(v) * v;

phi = DiracSpinorPlainWave(p, spinor_base)
phi = phi.';

%g5_phi = (g5 * phi).'; % Irrelevante

%cofs_g5 = DescomponSpinChiral(phi)
%phi_recompuesto = RecomponSpinChiral(cofs_g5)

% Ojo que con (1-g5) en vez de (Id-g5) no peta, pero da resultados malos.
Id = eye(4);
phi_L = .5*(Id-g5)*phi
phi_R = .5*(Id+g5)*phi


rho_L = phi_L' * phi_L
jx_L = (phi_L' * ax  * phi_L) / rho_L;
jy_L = (phi_L' * ay  * phi_L) / rho_L;
jz_L = (phi_L' * az  * phi_L) / rho_L;
j_L = [jx_L, jy_L, jz_L]
j_L_norm = sqrt(jx_L^2 + jy_L^2 + jz_L^2)

rho_R = phi_R' * phi_R
jx_R = (phi_R' * ax  * phi_R) / rho_R;
jy_R = (phi_R' * ay  * phi_R) / rho_R;
jz_R = (phi_R' * az  * phi_R) / rho_R;
j_R = [jx_R, jy_R, jz_R]
j_R_norm = sqrt(jx_R^2 + jy_R^2 + jz_R^2)

phi_i = gt * phi;
rho_i = phi_i' * phi_i
jx_i = (phi_i' * ax  * phi_i) / rho_i;
jy_i = (phi_i' * ay  * phi_i) / rho_i;
jz_i = (phi_i' * az  * phi_i) / rho_i;
j_i = [jx_i, jy_i, jz_i]

function cofs_g5 = DescomponSpinChiral(phi)
% Devuelve las 4 matrices gamma, ya que generalmente las necesito todas.

a1 = .5*(phi(1)+phi(3));
a2 = .5*(phi(1)-phi(3));

b1 = .5*(phi(2)+phi(4));
b2 = .5*(phi(2)-phi(4));

cofs_g5 = [a1, a2, b1, b2];

end

function f5_cofs = RecomponSpinChiral(cofs)

f5_cofs = [1,0,1,0] * cofs(1) + [1,0,-1,0] * cofs(2) + [0,1,0,1] * cofs(3) + + [0,1,0,-1] * cofs(4);

end

