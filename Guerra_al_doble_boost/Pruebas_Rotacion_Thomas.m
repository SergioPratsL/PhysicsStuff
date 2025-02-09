clear;

vA = [0.4, 0, 0];
vC = [0, 0.4, 0];

% Velocidad con la que A ve a B. El orden de los parámetros es al reves que
% los índices de v... así lo hice en su día
vAC = Vel_Addition_Law( vC, vA )

vCA = Vel_Addition_Law( vA, vC )

% La rotación se hace sobre el eje entre vA y vC
eje_rot = cross(vA, vC);
eje_rot = eje_rot / norm(eje_rot);
[e1, e2, e3] = ObtenBase(eje_rot);

angulo = acos(-dot(vAC, vCA) / (norm(vAC)*norm(vCA)));
rotMatrix = RotationMatrix(angulo);

% Cojo un vector espacial de B y lo envío a C pasando por A y aplicando para el
% boost menos la velocidad que ve C (vCA)

%test_vector_C = [1,0,0,0]
%test_vector_C = [0,1, 0, 1.71]
%test_vector_C = [0,0, -2, 1]
test_vector_C = [0,1.2, -2, 1.71]       % Funcionó!

test_vector_B = Boost(test_vector_C, -vC);

test_vector_A_dos_boosts_normal = Boost(test_vector_B, vA);          

test_vector_A_boost_directo_vCA = Boost(test_vector_C, vCA)

% Rotamos en destino el resultado de los dos boosts para compararlo con el
% resultado directo.
 test_vector_A_rotated = [test_vector_A_dos_boosts_normal(1), RotacionGeneral(test_vector_A_dos_boosts_normal(2:4), eje_rot, rotMatrix)]
 
 
% Ahora a ver cómo se hace si  usásemos -vAC en lugar de vCA
test_vector_C_rotated = [test_vector_C(1), RotacionGeneral(test_vector_C(2:4), eje_rot, rotMatrix)];

 
 test_vector_B_v2 = Boost(test_vector_C_rotated, -vC);

 test_vector_A_dos_boosts_rotated = Boost(test_vector_B_v2, vA)
 
 test_vector_A_boost_directo_vAC = Boost(test_vector_C, -vAC)
  
  
  % Ahora ver qué hay que hacer para que los boosts directos sean
  % iguales...
  
 test_vector_A_boost_directo_vCA_2 = Boost(test_vector_C_rotated, vCA);
 
 test_vector_A_boost_directo_vCA_2_rotated = [test_vector_A_boost_directo_vCA_2(1), RotacionGeneral(test_vector_A_boost_directo_vCA_2(2:4), eje_rot, rotMatrix')]
 

  % Arriba de eso no tocar!
 
 % Otras pruebas: veo que encaja con la velocidad propia, que es Sigma *
 % dv, y no hay rotación apreciable :).
 vLA = [0.8, 0, 0];
 vLB = [0.8, 0.001, 0];
 
 vect_t = [1, 0, 0, 0];
 vect_i_1 = Boost(vect_t, -vLA);
 vect_f_1 = Boost(vect_i_1, vLB)
 
 vect_i_2 = Boost(vect_t, -vLB);
 vect_f_2 = Boost(vect_i_2, vLA)
 
 
 vLA = [0.8, 0, 0];
 vLB = [0.001, 0, 0];
 
 % El factor de reducción es 1/Sigma^2
 vAB = Vel_Addition_Law(vLB, vLA);
 
 
 % Si sumas momento tenemos un factor 1/Sigma^3 pero  si sumas velocidades
 % El momento dp en origen se muevel Sigma*dp al hacer el boost. 
 vel = [0.8, 0, 0];
 p = fGamma(vel) * vel;
 p = p + [0.001, 0, 0];
 p_cuadr = norm(p)^2;
 Sigma = sqrt(1 + p_cuadr);
 vel2 = p / Sigma;
 % El factor de reducción es Sigma^3
 test = vel2 - vel
    
 