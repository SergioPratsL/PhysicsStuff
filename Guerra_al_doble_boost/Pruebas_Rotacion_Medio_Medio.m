% Este script tiene como objetivo ver que con diferentes velocidades la
% rotación medio medio lleva a matrices que son compatibles con un boost y
% para las cuales por fin vCA = -vAC completando lo que cerré en falso hace
% 15 años!

clear;

%vBA = [0.4, 0, 0];
%vBC = [0, 0.4, 0];
% vCA =     0.3837   -0.3837         0
% vAC =    -0.3837    0.3837         0

%vBA = [0.3, 0, 0];
%vBC = [0, 0.5, 0];
% vCA =   0.2801   -0.4889         0  
% vAC =  -0.2801    0.4889         0  

%vBA = [0.3, 0, -0.4];
%vBC = [0, 0.48, 0.1];
% vCA = 0.2737   -0.4366   -0.4558  
% vAC = -0.2737    0.4366    0.4558

%vBA = [0.79, -0.2, -0.4];
%vBC = [-0.15, 0.287, 0.3];
% vCA = 0.7360   -0.3383   -0.5094
% vAC = -0.7360    0.3383    0.5094

% Excelente, otro pequeño avance, parece que el giro medio medio es el daño
% de la rotación de Wigner y de la paradoja vCA + vAC...

% Más pruebas... parecido al MCU, ver la velocidad entre las dos posiciones
% cercanas, es decir, vAC
%vBA = [0.8, 0, 0];
%vBC = [0.8, 0.001, 0];
%vAC = - vCA = [6.9444e-07,    0.00166666         0]

% Ahora A pasa a ser B (el lab)
%vBA = [-0.8, 0, 0];
%vBC = [6.9444e-07, 0.00166666, 0];
% vAC = -vCA = 0.8000    0.0013         0

vBA = [-0.4, 0, 0];
vBC = [ -0.3837, 0.3837,  0];
% vAC = 0.0004    0.4159         0
% Los valores iniciales no vuelven, eso es inaceptable!!!



% Boostea de A a B
Boost_AB = Tensor_boosts(-vBA);
% Boostea de B a C
Boost_BC = Tensor_boosts(vBC);

% Boostea de C a B
Boost_CB = Tensor_boosts(-vBC);
% Boostea de B a A
Boost_BA = Tensor_boosts(vBA);

% Doble boost, él mismo no es un boost
Boost_AB_BC = Boost_AB * Boost_BC;
Boost_CB_BA = Boost_CB * Boost_BA;

vAC_AddLaw = Vel_Addition_Law( vBC, vBA )
vCA_AddLaw = Vel_Addition_Law( vBA, vBC )

eje_rot = cross(vBC, vBA);
angulo_medio = acos(-dot(vAC_AddLaw, vCA_AddLaw) / (norm(vAC_AddLaw)*norm(vCA_AddLaw))) / 2;
rotMatrix_media = RotationMatrixGeneral(eje_rot, angulo_medio);
rotMatrix4_media = RotMatrixTo4Rot(rotMatrix_media);

Boost_AB_BC_med_med = rotMatrix4_media * Boost_AB_BC * rotMatrix4_media
vCA = ObtenVelocidadDeMatrizBoost(Boost_AB_BC_med_med)

% Cambiar el orden de los vectores para que eje_rot apunte en sentido
% contrario y así la rotación sea al revés! (Se podría transponer pero me
% gusta más así)
eje_rot2 = cross(vBA, vBC);
rotMatrix_media2 = RotationMatrixGeneral(eje_rot2, angulo_medio);
rotMatrix4_media2 = RotMatrixTo4Rot(rotMatrix_media2);

Boost_CB_BA_med_med = rotMatrix4_media2 * Boost_CB_BA * rotMatrix4_media2
vAC = ObtenVelocidadDeMatrizBoost(Boost_CB_BA_med_med)

