
% Da el campo radiado por dos caminos, con la formula de Lienard-Wiechert con v <>0 y
% aplicando la misma formula para v=0 mas los boosts pertinentes
% Como siempre velocidad normalizada y las unidades brillan por su ausencia
% v: velocidad a la que se mueve la carga respecto al laboratorio.
% dv: aceleracion en el laboratorio
% R: distancia en el laboratorio

clear


% TRAS COMENTAR EN Boost_EM.m: 
% Eo = E_Sigma + BxV;
% Bo = B_Sigma - ExV;

%Prueba 1  fue OK (y vuelve a ser OK): OK con formula nueva
 % v = [0.4, 0, 0];
 % R = [1, 0, 0];
 % dv = [0, 1, 0];

%Prueba 2 (OK bajo ciertos criterios): (Es igual que la 1!!)
% v = [0.4, 0, 0];
% R = [-1, 0, 0];
% dv = [0, 1, 0];


% Prueba 3 (OK bajo mismos criterios que prueba 2): Ok con formula nueva
% v = [0.4, 0, 0];
% R = [0, 1, 0];
% dv = [-2, 2, 1];

% Prueba 4 (OK bajo los criterios anteriores): Ok con formula nueva
% v = [0.4, 0, 0];
% R = [1.5, 1, 0];
% dv = [0, -0.7, 1];

 
% Prueba 5 (OK): Ok con formula nueva
% v = [0.4, 0, 0];
% R = [1.2, -1, 0.9];
% dv = [-0.3, -0.7, 1.2];

 
% Prueba 6 (OK): Ok con formula nueva. Soy un maquina :D
 v = [-0.25, 0, 0.28];
 R = [1.2, -1, 0.9];
 dv = [-0.3, -0.7, 1.2];
 
 
[E_new, H_new] = CampoRadCargAcc_FormulaSergio( R, v, dv) 

[E_lab, H_lab] = CampoRadiadoCargaAcelerada( R, v, dv)


dist = norm(R);

% 4 vector de la luz
Rt = [dist, R(1), R(2), R(3)];

% Transformamos al sistema de la carga
Rt_carga = Boost( Rt, v );

R_carga = [Rt_carga(2), Rt_carga(3), Rt_carga(4)];

v0 = [0, 0, 0];

% Transformacion de la aceleracion propia
% norm_a = norm(a);
veloc = norm(v);
Sigma = fGamma(veloc);
%dv = a / Sigma;
e3 = v / veloc;
[e1, e2] = gPerpendicular(e3);

dv_new_base = CambioBase(e1, e2, e3, dv);

% La composicion de velocidades es clara (aunque a priori no lo parezca, 
% coordendas perpendiculares quedan divididas por Sigma, y la paralela por
% Sigma^2, aqui hacemos el camono a la inversa por lo que toca multiplicar

% Esta linea esta bien y debe mantenerse o si no las proporciones exactas
% se pierden, no caigamos en el panico de lo que ha fallado
dv_new_base_mod = [dv_new_base(1)*Sigma, dv_new_base(2)*Sigma, dv_new_base(3)*Sigma^2];
%dv_new_base_mod = dv_new_base;

% Recomponemos el vector
a = dv_new_base_mod(1)* e1 + dv_new_base_mod(2)* e2 + dv_new_base_mod(3)* e3;

% Esto es por la dilatancion del tiempo propio, antes lo tenia como 1/Sigma
% en el CampoRadiadoCargaAcelerada, pero tiene mas sentido aqui
a = a * Sigma;

% Campo que ve la carga en su sistema
[E_pre, H_pre] = CampoRadiadoCargaAcelerada( R_carga, v0, a);

% Ahora se ha de hacer un boost de vuelta al laboratorio,
% por tanto la velocidad sera -v!!
[Eo, Ho] = Boost_EM(E_pre, H_pre, -v)


% Obtenemos el efecto Doppler
% Dop = fDoppler(v, R);

% E_lab_bis = E_lab / Dop
% H_lab_bis = H_lab / Dop
