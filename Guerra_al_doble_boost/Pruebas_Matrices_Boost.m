clear;

vBA = [0.4, 0, 0];
vBC = [0, 0.4, 0];

% Velocidad con la que A ve a B. El orden de los parámetros es al reves que
% los índices de v... así lo hice en su día
vAC = Vel_Addition_Law( vBC, vBA )

vCA = Vel_Addition_Law( vBA, vBC )


% Boostea de A a B
Boost_AB = Tensor_boosts(-vBA);
% Boostea de B a C
Boost_BC = Tensor_boosts(vBC);

% Esta mierda no es simétrica, no es un boost.
Boost_AB_BC = Boost_AB * Boost_BC

% Boostea de C a A con -VCA.
Boost_CA_vCA = Tensor_boosts(vCA)      % ok
% Boostea de C a A con VAC.
Boost_CA_vAC = Tensor_boosts(-vAC)       % ok

% Boostea de A a C con -VAC.
Boost_AC_vAC = Tensor_boosts(vAC)      %
% Boostea de A a C con VCA.
Boost_AC_vCA = Tensor_boosts(-vCA);      %

eje_rot = cross(vBC, vBA);
angulo = acos(-dot(vAC, vCA) / (norm(vAC)*norm(vCA)));

rotMatrix = RotationMatrixGeneral(eje_rot, angulo);

rotMatrix4 = RotMatrixTo4Rot(rotMatrix);

Boost_AB_BC_rot = Boost_AB_BC * rotMatrix4

Dif = Boost_AB_BC_rot - Boost_AC_vAC;

angulo_medio = angulo / 2;
rotMatrix_media = RotationMatrixGeneral(eje_rot, angulo_medio);
rotMatrix4_media = RotMatrixTo4Rot(rotMatrix_media);
Boost_AB_BC_rot_new = rotMatrix4_media * Boost_AB_BC * rotMatrix4_media

v = ObtenVelocidadDeMatrizBoost(Boost_AB_BC_rot_new)

%v = [0.24, -0.53, 0.49];
%X = Tensor_boosts(-vAC)


