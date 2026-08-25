% addpath('C:\Fisica\Scripts matlab\HacerTuTeoria\8_origenes')

t = [1,0,0,0]

%v_BA = [0.4, 0, 0];
%v_BC = [0, 0.4, 0];

v_BA = [0.4, 0, 0];
v_BC = [-0.2, 0.5, 0.36];

B_BA = MatrizBoost(v_BA);
B_BAi = MatrizBoost(-v_BA);
B_BC = MatrizBoost(v_BC);
B_BCi = MatrizBoost(-v_BC);

v_tBA = t * B_BA;

B_CAdb = B_BA * B_BCi
%B_CAdbi = inv(B_CAdb);

v_CBA = t * B_CAdb;

% ACdb es transformación de C a A con doble boosts, es decir, C -> B -> A
B_ACdb = B_BC * B_BAi;
B_Acdbi = inv(B_ACdb);

% Da la identidad aunque son velocidades no inversas.
Expected_Id = B_CAdb * B_ACdb;



vect = t * B_CAdb;
veloc_CA = vect(2:4) / vect(1);

% Con velocidades inversas tiene que dar también la identidad.
v_CA = Vel_Addition_Law(v_BA, v_BC);


B_CAdir = MatrizBoost(v_CA); 
B_CAdiri = MatrizBoost(-v_CA); 

Otro_Expected_Id = B_CAdir * B_CAdiri;

v_AC = Vel_Addition_Law(v_BC, v_BA);
B_ACdir = MatrizBoost(v_AC); 

% La combinación del boost directo con el doble del sentido contrario
Not_Id = B_CAdir * B_ACdb;
Not_Id2 = B_ACdir * B_CAdb;


% Esta parte no sale bien.

rot = cross(v_BA, v_BC);
rot = rot / norm(rot);

[dir_rot, angulo] = Angulo_Wigner(v_BA, v_BC)

%[axis, angle] = Angulo_Wigner([0.4, 0, 0], [0, -0.4, 0]);

mat_rot = MatrizRotacion_4_4(dir_rot, angulo);
mat_rot_i = inv(mat_rot);


% check1 = inv(B_CAdir) - B_ACdir
% check2 = inv(mat_rot * B_CAdir) - B_ACdir
% % El peor de todos
% check3 = inv(mat_rot_i * B_CAdir) - B_ACdir
% 
% % Esta la peta:
% check41 = inv(mat_rot * B_CAdir * mat_rot_i) - B_ACdir
% check42 = inv(mat_rot_i * B_CAdir * mat_rot) - B_ACdir
% 
% check5 = inv(B_CAdir * mat_rot_i) - B_ACdir
% check6 = inv(B_CAdir * mat_rot) - B_ACdir
% 
% checkXX = inv(B_CAdb) - B_ACdir
% 
% check7 = inv(mat_rot * B_CAdb) - B_ACdir
% check8 = inv(mat_rot_i * B_CAdb) - B_ACdir
% check91 = inv(mat_rot * B_CAdb * mat_rot_i) - B_ACdir
% check92 = inv(mat_rot_i * B_CAdb * mat_rot) - B_ACdir
% check10 = inv(B_CAdb * mat_rot_i) - B_ACdir
% check11 = inv(B_CAdb * mat_rot) - B_ACdir

%check12 = B_CAdir * mat_rot - B_CAdb
%check13 = B_CAdir * mat_rot_i - B_ACdb

% Esto no sale ni para atrás.
%check14 = B_CAdir * mat_rot - mat_rot_i * inv(B_ACdir)
%check15 = B_CAdir * mat_rot_i - mat_rot * inv(B_ACdir)

%check16 = B_CAdir * mat_rot - mat_rot_i * B_ACdir
%check17 = B_CAdir * mat_rot_i - mat_rot * B_ACdir

B = sqrtm(B_CAdb * B_CAdb.');
% B es la inversa de B_ACdir, no es lo mismo que B_CAdir!
check20 = B - inv(B_ACdir)

R = inv(B) * B_CAdb;
% R es mat_rot, pero con velocidades altas empieza a fallar... pues lo que
% falla es el Angulo_Wigner, R es la verdadera rotación!
check21 = R - mat_rot


B_prim = sqrtm(B_CAdb.' * B_CAdb)
R_prim = B_CAdb * inv(B_prim) ;

% Todo Cuadra!
check_final = B * R - R_prim * B_prim
check_final2 = R - R_prim
check_final3 = B_prim - inv(R) * B * R

check_B_vs_B_dir1 = B - inv(B_ACdir)
check_B_vs_B_dir2 = B_prim - B_CAdir

% Todo ha quedado aclarado en el teorema de la descomposicion del doble
% boost
