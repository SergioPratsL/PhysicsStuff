u = [0.6, 0, 0];

x1 = [0, 1, 0];
v1 = [0.3, 0, 0];
x2 = [0,-1, 0];
v2 = [-0.3, 0, 0];

v1_lab = Vel_Addition_Law(v1, -u);
v2_lab = Vel_Addition_Law(v2, -u);

Sigma1 = fGamma(v1_lab);
Sigma2 = fGamma(v2_lab);

m_ang1 = Sigma1 * cross(x1, (v1_lab - u))
m_ang2 = Sigma2 * cross(x2, (v2_lab - u))

%m_ang1_no_internal = Sigma1 * cross(x1, (v1_lab))
%m_ang2_no_internal = Sigma2 * cross(x2, (v2_lab))

Sigma_lab = fGamma(u);

x3 = [-1, 0, 0];
v3 = [0, 0.3, 0];
v3_lab = Vel_Addition_Law(v3, -u);
Sigma3 = fGamma(v3_lab);
m_ang3 = Sigma3 * cross(x3, v3_lab) /  Sigma_lab

Sigma_ori = fGamma(v1);
m_ang_ori = Sigma_ori * cross(x1, v1)

u_new = [0, 0, 0.6];
x1_new = [0, 1, 0];
v1_new = [0.3, 0, 0];

v1_lab_new = Vel_Addition_Law(v1_new, -u_new);
Sigma1_new = fGamma(v1_lab_new);
m_ang1_new = Sigma1_new * cross(x1_new, (v1_lab_new - u_new))


% Cambiamos a u_new2, conservamos x1_new y v1_new
%u_new2 = [0, 0.4, sqrt(0.36-0.16)];
u_new2 = [0.4, 0, sqrt(0.36-0.16)];
v2_lab_new = Vel_Addition_Law(v1_new, -u_new2);
Sigma2_new = fGamma(v2_lab_new);

x1_new_contr = ContraeDireccion(x1_new, u_new2);

m_ang2_new = Sigma2_new * cross(x1_new_contr, (v2_lab_new - u_new2))

Sigma_z = 1/ fGamma([u_new2(1), u_new2(2),0]);
Ratio = norm(m_ang2_new) / norm(m_ang_ori);

angulo_rot = atan(m_ang2_new(1)/m_ang2_new(3));     % 0.1086

% Niquelado!
m_ang2_new_altern = ContraeDireccionPerpend(m_ang_ori, u_new2)