clear;

v_0 = [0,0,0];

% Prueba 1
%A = LeviCivita([1,2,2])
%B = LeviCivitaMatrix(3);
% C = LeviCivitaMatrix(4);
% D = sum(sum(sum(sum(abs(C)))));  % 24 valores no nulos de 256
%
% N = [0, 0, 0];
% J = [0, 0, .5];
% v = [0.4, 0, 0];
% 
% M_ori = TensorMomAng(N, J)
% M = Boost_Tensor_Momento_Angular(N, J, v)
% 
% P_ori = Momento_4(v_0);
% % Ojo a este signo, estamos evaluando alguien en el "0" desde alguien en
% % "v" por lo que la velocidad que ve es -v!
% %P = Momento_4(-v);
% 
% W_ori = Pseudovector_P_L(M_ori, P_ori)
% W = Pseudovector_P_L(M, P)
% %W = [0         0         0    0.5] Invariante OK


% Prueba 2
N = [0, 0, 0];
J = [0, 0, 0.5];
%v = [0, 0, 0.4];
v = [0, 0.4, 0];
%v = [0.4, 0.38, -0.27];

M_ori = TensorMomAng(N, J)
M = Boost_Tensor_Momento_Angular(N, J, v)

[N_boost, L_boost] = ExtraeEM_deTensor(M)

P_ori = Momento_4(v_0);

P = Momento_4(-v);

W_ori = Pseudovector_P_L(M_ori, P_ori)
W = Pseudovector_P_L(M, P)

% El vector P-L requiere boosts en sentido inverso al ser un vector axial.
W_boost = Boost(W_ori, -v)

%W = [0.2182         0         0    0.5455], Invariante ok
%W = [-0.1711   -0.0382   -0.0363    0.5258]


% Prueba 3, con momento de inercia arbitrario
% N = [1, 10, 100];
% %N = [0, 0, 0];
% J = [0, 0, 0.5];
% v = [0.4, 0.3, 0.5];
% 
% M_ori = TensorMomAng(N, J)
% M = Boost_Tensor_Momento_Angular(N, J, v)
% 
% % Adultermos el momento angular en otro SRI para comprobar que eso si hace
% % efecto. Hay que respetar la antisimetria.
% M_adulterado = M;
% M_adulterado(2,1) = M_adulterado(2,1) + 10;
% M_adulterado(1,2) = M_adulterado(1,2) - 10;
% 
% P_ori = Momento_4(v_0);
% P = Momento_4(-v);
% 
% W_ori = Pseudovector_P_L(M_ori, P_ori)
% W = Pseudovector_P_L(M, P)
% W_adulterado = Pseudovector_P_L(M_adulterado, P)
% inv_W_ori = - Distancia4Vector(W_ori)
% inv_W = - Distancia4Vector(W)
% inv_W_adulterado = - Distancia4Vector(W_adulterado)
% % W = 0.3536    0.0828    0.0621    0.6036  con N
% % W = 0.3536    0.0828    0.0621    0.6036 sin N (no hace nada!)
% % W_adulterado = 0.3536    0.0828    7.1332   -3.6391
% % Cuando se está fuera del SRI origen el momento adulterado sí que afecta
% % al vector P-L.
% 
% % Si hay movimiento, las componentes de N perpendiculares a la velocidad
% % dan aportaciones al vector PL (a sus componentes espaciales).





% Prueba 4, pensando que la partícula tiebe OAM.
% Supongamos que L=1 y S=0, voy a necesitar librarme del spin por un
% momento, paree que da igual que la partícula se mueva, si quiero
% sobrevivir tengo que coger el P del centro de masa, es decir el origen P
% = [1,0,0,0] aunque la partícula se mueva, lo contrario es la muerte!

% v_ori = [0, 0.4, 0];
% %P_ori = Momento_4(v_ori);
% P_ori = Momento_4(v_0);       % v_ori debe ser ignorada!
% 
% % En verdad N origen es irrelevante :).
% N = [0, 0, 0];
% J = [0, 0, 1];
% v_boost = [0, 0, 0.4];     
% 
% M_ori = TensorMomAng(N, J);
% M = Boost_Tensor_Momento_Angular(N, J, v_boost);
% 
% % v_new_sri = Vel_Addition_Law(v_ori, v_boost);   % deshacerse de v_ori!
% v_new_sri = Vel_Addition_Law(v_0, v_boost);
% 
% P = Momento_4(v_new_sri);
%  
% W_ori = Pseudovector_P_L(M_ori, P_ori)
% W = Pseudovector_P_L(M, P)
% % Un poco redundante el cambio de nomenclatura de esta prueba, pero es que
% % me dio problemas decidir qué debia ser P_ori y demás.
% % Dió bien, lo esperado.



% Ver si los boosts pueden transformar directmaente en vector P-L
% N = [0, 0, 0];
% %J = [10, 100, 1000];
% J = [2, -5, 9];
% v = [0.2, -0.3, 0.4];
% %v = [0.4, 0, -0.1];
% 
% M_ori = TensorMomAng(N, J)
% M = Boost_Tensor_Momento_Angular(N, J, v)
% 
% P_ori = Momento_4(v_0);
% % Ojo a este signo, estamos evaluando alguien en el "0" desde alguien en
% % "v" por lo que la velocidad que ve es -v!
% P = Momento_4(-v);
% 
% W_ori = Pseudovector_P_L(M_ori, P_ori)
% W = Pseudovector_P_L(M, P)
% 
% W_boost = Boost(W_ori, -v)
% El boost tuvo todas las componentes iguales pero es extraño que tengo que
% usar -v en el boost cuando le correspondería usaar v. P es lo que ve el
% nuevo SRI respecto al viejo, pero el boost es del SRI origen al nuevo

