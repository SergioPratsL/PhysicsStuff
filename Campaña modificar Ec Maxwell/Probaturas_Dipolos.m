nivel_probaturas = 9.7


% Probar primero de todo que para un dipolo eléctrico en movimiento pequeño
% pero no infinitesimal, el dipolo se orientará respecto a la fuerza de
% Lorentz, y eso es equivalente a sacar el momento magnético y  aplicarle
% a tal dipolo magnético el magnetismo... va a ser que no :P
if nivel_probaturas == 1
   % Ante los fracasos le doy la vuelta y parto del SRI origen
   % u es el boost desde el SRI de la partícula al laboratorio  y 
   % E y B son los campos vistos en el SRI de la partícula
   u = [0.4, 0, 0];
   E = [1, 0.3, 0];
   B = [0, 0, 0];
   
   
   % La longitud del dipolo en el sistema origen es 1:
   l_origin = 1;
   
   Sigma = 1 / (1 - norm(u)^2)^(1/2);
   
   T_boost = Tensor_boosts(u);
   
   u_norm = u / norm(u);
   
   % Para determinar orientacion dipolo en el SRI original debo sabe 
   % hacia donde tira el campo E en dicho sistema
   [E_lab, B_lab] = Boost_EM(E, B, u);
   
   % Probatura sobre el vector origen
   cos_rot = 1.09 / (norm([1, 0.3*Sigma])*norm([1, 0.3/Sigma]));
   rot = -acos(cos_rot);        % Probar con los dos signos
   E_aux = [0,0, 0];
   E_aux(1) = cos(rot) * E(1) - sin(rot) * E(2);
   E_aux(2) = cos(rot) * E(2) + sin(rot) * E(1);
   
   % Asumo que el "dipolo absoluto" es 1, y debe medirse en el SRI origen
   Dip_E_ori = E / norm(E);
   Dip_E_ori = E_aux / norm(E_aux);     % Probatura sobre el vector origen
   % Es un barra eléctrica, en su sistema origina no tiene dip magnetico
   Dip_M_ori = [0, 0, 0];
     
   % Segun Boost_EM, veamos que direccion debería tener el dipolo ahora.
   %[Dip_E_boost, Dip_M_boost] = Boost_EM( Dip_E_ori, Dip_M_ori, -v_dip)
   % Forma correcta de hacer boosts de matrices.
   T_dip_ori = TensorEM(Dip_E_ori, Dip_M_ori);
   
   Invariante_dip_ori = norm(Dip_E_ori)^2 - norm(Dip_M_ori)^2;
   
   %T_dip_lab_boost = (T_boost * ((T_boost * T_dip_ori)'))'
   T_dip_lab_boost = (T_boost * T_dip_ori) * T_boost'      % Equivalente
   
   % Probatura de la rotación
   %T_aux = (T_boost * T_dip_ori);
   
   %T_rot = zeros(4);
   %T_rot(1,1) = 1;
   %T_rot(4,4) = 1;
   %T_rot(2,2) = cos(rot);
   %T_rot(3,3) = cos(rot);
   %T_rot(3,2) = sin(rot)
   %T_rot(2,3) = -sin(rot);
   
   %T_aux = T_aux * T_rot'; 
   %T_dip_lab_boost = T_aux  * T_boost'
   
   %T_dip_ori_aux = T_dip_ori * T_rot;
   %T_dip_lab_boost = (T_boost * T_dip_ori_aux) * T_boost';      % Equivalente
   
   [Dip_E_lab_boost, Dip_M_lab_boost] = SacaCamposTensorEM(T_dip_lab_boost);
   Invariante_dip_lab_boost = norm(Dip_E_lab_boost)^2 - norm(Dip_M_lab_boost)^2
   
   % Veamos por otro lado cómo es la fuerza de Lorent que indica hacia
   % donde ha de quedar orientado el dipolo
   F = FLorentz(-u, E_lab, B_lab);
   
   % Probatura
   p_ori = [0, E];
   p_lab_boost = Boost(p_ori, u);
   
   Dir_F = F / norm(F);

   comp_dir_veloc = dot( Dir_F, u_norm);
   
   % Calcularlo de otra forma
   %l = l_origin;
   %l = l * sqrt( (comp_dir_veloc / Sigma)^2 + (1 - comp_dir_veloc^2));
   val_no_contract = sqrt( (Sigma * dot(F, u_norm))^2 + (norm(F)^2 - dot(F, u_norm)^2) );
   % Obtenemos que proporción del total se contrae la barra, reduciendo el
   % momento
   coef_contract = norm(F) / val_no_contract;
   l = l_origin * coef_contract;
   
   
   Dip_E_lab = Dir_F * l;
   
   Dip_m_lab = cross( Dir_F, -u) * l;
   
   T_dip_lab_direct = TensorEM(Dip_E_lab, Dip_m_lab);
   % Factor de aumento de densidad debido a la contracción de una dimensión
   T_dip_lab_direct = T_dip_lab_direct * Sigma
   

   [Dip_E_lab_direct, Dip_M_lab_direct] = SacaCamposTensorEM(T_dip_lab_direct);
   Invariante_dip_lab_direct = norm(Dip_E_lab_direct)^2 - norm(Dip_M_lab_direct)^2;
   
   Dip_E_lab_direct = Dip_E_lab_direct;
   Dip_E_lab_boost = Dip_E_lab_boost;
   % Como esperaba da lo mismo que el boost.
   Dip_E_test = Sigma * (Dip_E_ori + cross(Dip_M_ori, u)) - (Sigma - 1) * dot(u_norm, Dip_E_ori) * u_norm;
   

% La forma en que transforma la orientación del dipolo es el resultado de
% hacer un boost tanto a los vectores del dipolo como al vector temporal
% al cual se le debe hacer un boost.
elseif nivel_probaturas == 2
   v_boost = [0.4, 0, 0];
   E = [1, 0.3, 0];
   B = [0, 0, 0];
   %B = [0, -0.25, 0];
   
   Sigma = 1 / (1 - norm(v_boost)^2)^(1/2);
   v_boost_norm = v_boost / norm(v_boost);
   l_origin = 1;
   
   % Asumo que el "dipolo absoluto" es 1, y debe medirse en el SRI origen
   Dip_E_ori = E / norm(E);
   % Es un barra eléctrica, en su sistema origina no tiene dip magnetico
   Dip_M_ori = [0, 0, 0];
     
   T_dip_ori = TensorEM(Dip_E_ori, Dip_M_ori)
   
   
   u_ori = [1, 0, 0, 0];
   T_EM_ori = TensorEM(E, B);
   
   % Con T_EM_ori' sale la fuerza con signo inverso!
   F_ori = u_ori * T_EM_ori';
   dir_dip_ori = F_ori / norm(F_ori);
   
   T_dip_ori_new_deal = u_ori' * dir_dip_ori - (u_ori' * dir_dip_ori)';
   
   [E_lab, B_lab] = Boost_EM(E, B, v_boost);
   % Probatura E en vez de F: orientar el dipolo en función del campo
   % eléctrico y no del campo total
   %F_lab = FLorentz(-v_boost, E_lab, B_lab);

   val_no_contract = sqrt( (Sigma * dot(E_lab, v_boost_norm))^2 + (norm(E_lab)^2 - dot(E_lab, v_boost_norm)^2) );
   coef_contract = norm(E_lab) / val_no_contract;
   
   
   l = l_origin * coef_contract;
   Dir_E = E_lab / norm(E_lab);
   Dip_E_lab = Dir_E * l;
   Dip_m_lab = cross( Dir_E, -v_boost) * l;
   
   T_dip_lab_direct = TensorEM(Dip_E_lab, Dip_m_lab);
   % Factor de aumento de densidad debido a la contracción de una dimensión
   T_dip_lab_direct = T_dip_lab_direct * Sigma
   
   u_lab = Boost(u_ori, v_boost);
   % las componentes 2 y 3 muestran misma relaciób que F_lab
   dir_dip_lab = Boost(dir_dip_ori, v_boost);
   
   % Lleva exactamente al mismo valor!
   % u_lab es un boost a un vector temporal.
   T_dip_lab_new_deal = (u_lab' * dir_dip_lab)' - (u_lab' * dir_dip_lab);
   

   T_boost = Tensor_boosts(v_boost);
   T_dip_lab_boost = (T_boost * T_dip_ori) * T_boost'  
   
   
   [Dip_E_lab_new_deal, Dip_M_lab_new_deal] = SacaCamposTensorEM(T_dip_lab_new_deal);
   Invariante_dip_lab_new_deal = norm(Dip_E_lab_new_deal)^2 - norm(Dip_M_lab_new_deal)^2;

   
% Comparar cómo se transforma la dirección del bastón con cómo se
% transforma el momento dipolar. El baston se transforma siguiendo la norma
% de la contracción de la longitud paralela a la velocidad
elseif nivel_probaturas == 3   
   % Boost siempre en dirección X para simplificar algunas cosas!
   v_boost = [0.4, 0, 0];       
   % Test 1 OK.
   E = [1, 0.3, 0];
   B = [0, 0, 0];
   % Test 2 OK.
   E = [0.4, -0.3, 0.5];
   B = [0, 0, 0];
   % Test 3 OK
   E = [1.1, -0.3, 0.5];
   B = [-0.2, 0.7, 0];
   % Test 4 OK
   E = [1.1, 0.08, 0.5];
   B = [-0.2, 0.7, 0.6];
   
   Sigma = 1 / (1 - norm(v_boost)^2)^(1/2);
   v_boost_norm = v_boost / norm(v_boost);
   l_origin = 1;
   
   Dip_E_ori = E / norm(E);
   Dip_M_ori = [0, 0, 0];
   
   Barra_lab_boost = [Dip_E_ori(1) / Sigma, Dip_E_ori(2), Dip_E_ori(3)];
   Barra_lab_boost_norm = Barra_lab_boost / norm(Barra_lab_boost)
   
   T_dip_ori = TensorEM(Dip_E_ori, Dip_M_ori);
   
   T_boost = Tensor_boosts(v_boost);
   
   T_dip_lab_boost = (T_boost * T_dip_ori) * T_boost' ; 
   
   Barra_lab_transf = T_dip_lab_boost(2:4, 1)';
   Barra_lab_transf_norm = Barra_lab_transf / norm(Barra_lab_transf)
   
   
% Cancelado
elseif nivel_probaturas == 4
    v_boost = [0.4, 0, 0];       
    vel = norm(v_boost);
    Sigma = 1 / (1 - norm(v_boost)^2)^(1/2);
    
    % Misma Y, diferente X, formar un triangulo
    A = [Sigma, 0.3];
    B = [1/Sigma, 0.3];
    C = Sigma * [vel^2, 0, 0];       % ???
    % Cuando más vertical sea la dirección del bastón, más pequeño será C
    % en relción con la longitud total del bastón.
   
% Ver si rotando la dirección de la fuerza de la misma forma que se rota
% "el acelerin", podemos conseguir algo.
elseif nivel_probaturas == 5    
   v_boost = [0.4, 0, 0];       
   vel = norm(v_boost);
   Sigma = 1 / (1 - norm(v_boost)^2)^(1/2);
    
   E = [1, 0.3, 0];
   B = [0, 0, 0];
   
   F_ori = E;
   F_ori_norm = F_ori / norm(F_ori)
    
   Sigma = 1 / (1 - norm(v_boost)^2)^(1/2);
   v_boost_norm = v_boost / norm(v_boost);
   l_origin = 1;
   
   Dip_E_ori = E / norm(E);
   Dip_M_ori = [0, 0, 0];    

   T_dip_ori = TensorEM(Dip_E_ori, Dip_M_ori);
   T_boost = Tensor_boosts(v_boost);   

   T_dip_lab_boost = (T_boost * T_dip_ori) * T_boost' ; 
   
   Barra_lab_transf = T_dip_lab_boost(2:4, 1)';
   Barra_lab_transf_norm = Barra_lab_transf / norm(Barra_lab_transf);   
   
   [E_lab, B_lab] = Boost_EM(E, B, v_boost);
   F_lab = FLorentz(-v_boost, E_lab, B_lab);
   F_lab_norm = F_lab / norm(F_lab)

   % Variacion de la velocidad vista desde el origen para calcular el giro
   dv_ori = [F_ori(1) / Sigma^2, F_ori(2), F_ori(3)];
   dv_ori = dv_ori / norm(dv_ori);
   
   dv_lab = 1/Sigma * [F_lab(1)/Sigma, F_lab(2), F_lab(3)];
   dv_lab = dv_lab / norm(dv_lab);
   
   cos_fase = dot(dv_ori, dv_lab);
   fase = acos(cos_fase);
   
   % Lo que hay rotar es la dirección de la barra en el SRI del 
   % laboratorio, no la fuerza original... para ver si coincide con la
   % fuerza
   F_ori_rot_A = Rota_puto_vector(Barra_lab_transf_norm, fase)
   F_ori_rot_B = Rota_puto_vector(Barra_lab_transf_norm, -fase)
   

% Demostrar que en Lab la en la dirección del baston, la dirección del 
% campo que viene del otro origen es la misma que la fuerza EM que se ve
% en el laboratorio
elseif nivel_probaturas == 6    
   v_boost = [0.4, 0, 0];       
   vel = norm(v_boost);
   Sigma = 1 / (1 - norm(v_boost)^2)^(1/2);
   % Test 1. OK  
   E = [1, 0.3, 0];
   B = [0, 0, 0];
   % Test 2. OK
   E = [0.4, -0.3, 0.5];
   B = [0, 0, 0];
   % Test 3. OK
   E = [1.1, -0.3, 0.5];
   B = [-0.2, 0.7, 0];
   % Test 4 
   E = [1.1, 0.08, 0.5];
   B = [-0.2, 0.7, 0.6];   
   
   F_ori = E;
   F_ori_norm = F_ori / norm(F_ori)
    
   v_boost_norm = v_boost / norm(v_boost);
   l_origin = 1;
   
   Dip_E_ori = E / norm(E);
   Dip_M_ori = [0, 0, 0];    

   T_dip_ori = TensorEM(Dip_E_ori, Dip_M_ori);
   T_boost = Tensor_boosts(v_boost);   

   T_dip_lab_boost = (T_boost * T_dip_ori) * T_boost' ; 
   
   Barra_lab_transf = T_dip_lab_boost(2:4, 1)';
   Rod = Barra_lab_transf / norm(Barra_lab_transf)
   
   % Sacar la distancia retardada
   cof_c = (Rod(1)^2 - vel^2*(Rod(2)^2 + Rod(3)^2));
   polin = [(1-vel^2), -2*Rod(1), cof_c];
   sols = roots(polin) 
   
   D1 = [sols(1), Rod(2), Rod(3)];
   dist_D1 = norm(D1);
   dif_1 = abs( D1(1) - (Rod(1) - v_boost(1) * dist_D1) );
   
   D2 = [sols(2), Rod(2), Rod(3)];
   dist_D2 = norm(D2);
   dif_2 = abs( D2(1) - (Rod(1) - v_boost(1) * dist_D2) );
    
   if dif_1 < dif_2
       R_ret = D1;
   else
       R_ret = D2;
   end
   
   % Calcula el campo EM causado por la fuerza en R_ret
   [Eo, Bo] = CampoInducido_sin_unidades(R_ret, -v_boost);
   
   Fo = FLorentz(-v_boost, Eo, Bo);
   Fo_norm = Fo / norm(Fo)
   
   [E_lab, B_lab] = Boost_EM(E, B, v_boost);
   F_lab = FLorentz(-v_boost, E_lab, B_lab);
   F_lab_norm = F_lab / norm(F_lab)
   
   
% Covarianza de la fuerza hecha sobre un dipolo tipo barra infinitesimal
% Entran en juego las derivadas del campo EM... empezaré jugando con las
% matrices de derivadas
% 13/05/2018. Demostrado que la fuerza en un mismo instante en el SRI
% propio y en el SRI externo no puede ser covariante cuando hay derivadas 
% temporales del campo, miraré de ver que la fuerza en diferentes instantes 
% del sistema externo (léase "lab") es la covariante.
elseif nivel_probaturas == 7       
   v_boost = [0.4, 0, 0];       
   vel = norm(v_boost);
   Sigma = 1 / (1 - norm(v_boost)^2)^(1/2);
   
   % Prueba 1. 
   %E = [1, 0, 0];
   %B = [0, 0, 0];
   %Cofs = [0.5, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; 
   
   % Prueba 2. OK
   %E = [1, 0.3, 0];
   %B = [0, 0, 0];
   %Cofs = [0.5, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; 
   
   % Prueba 3. OK
   %E = [1, 0.3, 0];
   %B = [0.1, -0.25, 0.18];
   %Cofs = [0.5, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
   
   % Prueba 4. KO
   %E = [1, 0, 0];
   %B = [0, 0, 0];
   %Cofs = [0.5, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2, -0.15, 0, 0]; 
   
   % Prueba 4B. también falla
   E = [1, 0.3, 0];
   B = [0.1, -0.25, 0.18];
   Cofs = [0.5, -1, 0.12, -0.08, -0.14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];       
   Cofs = [0.5, -1, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];       
   
   % [dxEx, dxEy, dtEx, dtEy, dtEz, dxEy, dxEz, dyEx
   % dyEz, dzEx, dzEy, dxBx, dyBy, dxBy, dyBz, dzBx]
   
   % Ojo! F_ori es la fuerza que orienta el dipolo, no hace fuerza efectiva
   F_ori = E;
   F_ori_norm = F_ori / norm(F_ori)
   
   l_origin = 1;
   Dip_E_ori = F_ori_norm;
   Dip_M_ori = [0, 0, 0];       
   
   l_origin = 1;
   Dip_E_ori = E / norm(E);
   Dip_M_ori = [0, 0, 0];    
   % 4 vector espacial, primera componente siempre cero!
   Barra_ori = [0, Dip_E_ori];
   
   % Transformar el momento    
   [E_uv, B_uv] = fCrea_Matriz_Derivadas(Cofs);
   F_deriv_ori = Barra_ori * E_uv;
   Four_mom_ori = [0, F_deriv_ori(2:4)]
   
   Four_mom_lab_boost = Boost(Four_mom_ori, v_boost)

   
   T_dip_ori = TensorEM(Dip_E_ori, Dip_M_ori);
   T_boost = Tensor_boosts(v_boost);   
   T_dip_lab_boost = (T_boost * T_dip_ori) * T_boost' ; 
   
   
   [Dip_E_lab, Dip_M_lab] = ObtenCampoEMTensor(T_dip_lab_boost);

   % MAAAALLLL
   %E_uv_lab = (T_boost * E_uv) * T_boost' ; 
   %B_uv_lab = (T_boost * B_uv) * T_boost' ; 
   [E_uv_lab, B_uv_lab] = Boost_Derivadas_EM(E_uv, B_uv, v_boost);
   
   Barra_lab = [0, Dip_E_lab];
   
   F_deriv_elec_lab = Barra_lab * E_uv_lab;
   F_deriv_elec_lab = F_deriv_elec_lab(2:4);
   B_dir_barra = Barra_lab * B_uv_lab;
   F_deriv_mag_lab = cross(-v_boost, B_dir_barra(2:4));
   
   F_deriv_lab = F_deriv_elec_lab + F_deriv_mag_lab;
   
   % Obtener el 4 momento. Dilatación temporal y contracción de longitudes
   % se compensan!!
   Four_mom_lab_direct = [dot(-v_boost, F_deriv_lab), F_deriv_lab] 
   
% 13/05/2018. Demostrado que la fuerza en un mismo instante en el SRI
% propio y en el SRI externo no puede ser covariante cuando hay derivadas 
% temporales del campo, miraré de ver que la fuerza en diferentes instantes 
% del sistema externo (léase "lab") es la covariante.
% El sufijo "lab" pasa a ser "ext" de externo, la razón es que lo lógico
% sería tener el dipolo en el laboratorio (ori) y luego evaluarlo desde
% sistemas a gran velocidad que realmente no existen (pero podrían existir)
elseif nivel_probaturas == 8       
   v_boost = [0.4, 0, 0];       
   vel = norm(v_boost);
   Sigma = 1 / (1 - norm(v_boost)^2)^(1/2);
   
   % Prueba 1. OK
   %E = [1, 0, 0];
   %B = [0, 0, 0];
   %Cofs = [0.5, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; 
   
   % Prueba 2. OK
   %E = [1, 0.3, 0];
   %B = [0, 0, 0];
   %Cofs = [0.5, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; 
   
   % Prueba 3. OK
   %E = [1, 0.3, 0];
   %B = [0.1, -0.25, 0.18];
   %Cofs = [0.5, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
   
   % Prueba 4. OK
   %E = [1, 0.3, 0];
   %B = [0, 0, 0];
   %Cofs = [0.5, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2, -0.15, 0, 0]; 
   
   % Prueba 5. OK
   %E = [1, 0.3, 0];
   %B = [0.1, -0.25, 0.18];
   %Cofs = [0.5, -1, 0.12, -0.08, -0.14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];       
   
   % Prueba 6. OK
   %E = [1, 0.3, -0.25];
   %B = [0.1, -0.25, 0.18];
   %Cofs = [0.5, -1, 0.12, -0.08, -0.14, 0.2, -0.44, 0, 0, 0, 0, 0, 0, 0, 0, 0];       
   
   % Prueba 7. OK
   %E = [1, 0.5, -0.25];
   %B = [0.1, 0, 0.18];
   %Cofs = [0.5, -1, 0.12, -0.08, -0.14, 0.2, -0.44, 0.1, -0.32, 0.23, 0, 0, 0, 0, 0, 0];       
   
   % Prueba 8. Victoria!!!!!
   E = [1, 1.5, -0.25];
   B = [1, -0.3, 0.18];
   Cofs = [0.5, -1, 0.12, -0.08, -0.14, 0.2, -0.44, 0.1, -0.32, 0.23, -0.34, 1, 0.3, -0.7, 0.28, 0.45];       
      
   
   % [dxEx, dxEy, dtEx, dtEy, dtEz, dxEy, dxEz, dyEx
   % dyEz, dzEx, dzEy, dxBx, dyBy, dxBy, dyBz, dzBx]
   
   % Ojo! F_ori es la fuerza que orienta el dipolo, no hace fuerza efectiva
   F_ori = E;
   F_ori_norm = F_ori / norm(F_ori)
   
   l_origin = 1;   
   Dip_E_ori = E / norm(E);
   Dip_M_ori = [0, 0, 0];    
   % 4 vector espacial, primera componente siempre cero!
   Barra_ori = [0, Dip_E_ori];
   
   % Transformar el momento    
   [E_uv, B_uv] = fCrea_Matriz_Derivadas(Cofs);
   F_deriv_ori = Barra_ori * E_uv;
   Four_mom_ori = [0, F_deriv_ori(2:4)]
   
   Four_mom_ext_boost = Boost(Four_mom_ori, v_boost)
    
   T_boost = Tensor_boosts(v_boost);   
   
   %T_dip_ori = TensorEM(Dip_E_ori, Dip_M_ori);  
   %T_dip_lab_boost = (T_boost * T_dip_ori) * T_boost' ;   
   %[Dip_E_lab, Dip_M_lab] = ObtenCampoEMTensor(T_dip_lab_boost);

   [E_uv_ext, B_uv_ext] = Boost_Derivadas_EM(E_uv, B_uv, v_boost);
   
   % Esto ya no puede hacerse así, covarianza imposible de esta manera
   %Barra_lab = [0, Dip_E_lab];
   Barra_ext = Boost(Barra_ori, v_boost);
   
   F_deriv_elec_ext = Barra_ext * E_uv_ext;
   F_deriv_elec_ext = F_deriv_elec_ext(2:4);
   B_dir_barra = Barra_ext * B_uv_ext;
   F_deriv_mag_ext = cross(-v_boost, B_dir_barra(2:4));
   
   F_deriv_ext = F_deriv_elec_ext + F_deriv_mag_ext;
   
   % Obtener el 4 momento aplicando la dilatación temporal.
   Four_mom_ext_direct = [dot(-v_boost, F_deriv_ext), F_deriv_ext] * Sigma
   
   T_boost = Tensor_boosts(v_boost);   
   T_dip_ori = TensorEM(Dip_E_ori, Dip_M_ori);  
   T_dip_lab_boost = (T_boost * T_dip_ori) * T_boost' ;   
   T_dip_lab_boost = T_dip_lab_boost';
   [Dip_E_lab, Dip_M_lab] = ObtenCampoEMTensor(T_dip_lab_boost);
   asdade = 1;   % Párame el debug!

% Comprobar que la fuerza total sobre un dipolo equilibrado sigue siendo
% cero al cambiar de SRI.
elseif nivel_probaturas == 9       
   v_boost = [0.4, 0, 0];      
   vel = norm(v_boost);
   v_norm = v_boost / vel;
   Sigma = 1 / (1 - norm(v_boost)^2)^(1/2);
   
   % Prueba 1. OK
   %E = [0, 0, 0];
   %B = [0, 1, 0];

   % Prueba 2. OK
   %E = [0, 0, 0];
   %B = [-3, 1, 0];   
   
   % Prueba 3. OK
   %E = [0, -1.5, 0];
   %B = [-3, 1, -0.8];   
   
   % Prueba 4. OK
   E = [2, -1.5, 4.8];
   B = [5, -3, -1.8];      
   
   % Prueba 5. OK. Victoria... Pírrica.
   E = [6, -2.5, 4.8];
   B = [-2, 0, -2.8];         
   
   % Correccion a la fuerza eléctrica que obviamente está descompensada 
   % orientes hacia donde orientes el dipolo :).
   F_E_ori = 4 * E;
   F4_E_lab = Boost( [0, F_E_ori], v_boost);
   F_E_lab = F4_E_lab(2:4);
   
   momento_dipolo_mag_ori = 1;
   Dip_E_ori = [0, 0, 0];
   Dip_M_ori = B / norm(B) * momento_dipolo_mag_ori;       
   Dip_M_ori_norm = Dip_M_ori / norm(Dip_M_ori);
   
   %F_ori = 
   
   % Primero determinar la base para construir la espira cuadrada, hay que
   % cambiar de base. Espira en e2 (parcialmente paralela a la velocidad)
   % y en e3 que es perpendicular.
   [e1, e2, e3] = ObtenBaseEnh(Dip_M_ori/norm(Dip_M_ori), -v_norm);   
   
  
   % Obtener las velocidades en cada uno de los lados vistos desde el SRI
   % destino
   % lados [parcialmente] paralelos al boost al ver desde el
   % laboratorio (i.e. sistema externo)
   v1 = Vel_Addition_Law( e2, v_boost);
   v2 = Vel_Addition_Law( -e2, v_boost);
   % lados perpendiculares al boost
   v3 = Vel_Addition_Law( e3, v_boost);
   v4 = Vel_Addition_Law( -e3, v_boost);

   % Alargamientos de cada lado: creo que finalmente he dado con la
   % solución, si no es así apago y me voy
   long_x = abs(e2(1)) / Sigma;
   if long_x == 0
       l1 = Sigma;
       l2 = Sigma;
   else
       l1 = long_x / abs(v1(1) + v_boost(1));   % El lado se escapa a -v_boost
       l2 = long_x / abs(v2(1) + v_boost(1));
   end
   
   l3 = Sigma;
   l4 = Sigma;
   
   
   Dopler1 = Sigma * (1 - dot(e2, v_boost));
   Dopler2 = Sigma * (1 - dot(-e2, v_boost));
   Dopler3 = Sigma * (1 - dot(e3, v_boost));
   Dopler4 = Sigma * (1 - dot(-e3, v_boost));
   
   
   Dopler1_aux = 1/Sigma / (1 - dot(v1, -v_boost));
   Dopler2_aux = 1/Sigma / (1 - dot(v2, -v_boost));
   Dopler3_aux = 1/Sigma / (1 - dot(v3, -v_boost));
   Dopler4_aux = 1/Sigma / (1 - dot(v4, -v_boost));   
   
   % -Sigma * v_boost(1) * 4; es la distancia que la partícula se ha movido
   % en el sistema del laboratorio, por ello debe ser compensado para
   % verificar que cuadra
   Recorrido_Total_X = l1 * v1(1)+ l2 * v2(1) + l3*v3(1) + l4*v4(1) + Sigma * v_boost(1) * 4
   Recorrido_Total_Y = l1 * v1(2)+ l2 * v2(2) + l3*v3(2) + l4*v4(2)
   Recorrido_Total_Z = l1 * v1(3)+ l2 * v2(3) + l3*v3(3) + l4*v4(3)
   
   [E_lab, B_lab] = Boost_EM(E, B, v_boost);
   
   % Las l's y los doplers son lo mismo!
   F_mag1 = cross(v1, B_lab) * l1;
   F_mag2 = cross(v2, B_lab) * l2;
   F_mag3 = cross(v3, B_lab) * l3;
   F_mag4 = cross(v4, B_lab) * l4;
   
   F_mag = F_mag1 + F_mag2 + F_mag3 + F_mag4
   
   F_tot1 = F_mag1 + E_lab * l1;
   F_tot2 = F_mag2 + E_lab * l2;
   F_tot3 = F_mag3 + E_lab * l3;
   F_tot4 = F_mag4 + E_lab * l4;
   
   F_tot = F_tot1 + F_tot2 + F_tot3 + F_tot4 - F_E_lab
   
   %T_boost = Tensor_boosts(v_boost);   
   %T_dip_ori = TensorEM(Dip_E_ori, Dip_M_ori);  
   %T_dip_lab_boost = (T_boost * T_dip_ori) * T_boost' ;   
   %T_dip_lab_boost = T_dip_lab_boost';
   %[Dip_E_lab, Dip_M_lab] = ObtenCampoEMTensor(T_dip_lab_boost);

% Prueba de soporte al nivel 9, a pesar de la victoria se constata que una 
% corriente sin carga genera "fuerzas de mierda" es decir, fuerzas no 
% covariantes y por tanto reprobables y reprobadas desde el punto de vista
% física, es más con el estad ode arte actual se merecen el brutal
% calificativo de fuerzas "no-físicas"
elseif nivel_probaturas == 9.5   
   v_boost = [0.4, 0, 0];      
   vel = norm(v_boost);
   v_norm = v_boost / vel;
   Sigma = 1 / (1 - norm(v_boost)^2)^(1/2);
   
   % Prueba 1. OK 
   %E = [0, 0, 0];
   %B = [0, 1, 0];

   % Prueba 2. OK
   %E = [0, 0, 1];
   %B = [0, -3, 0];   
    
   % Prueba 3. OK
   %E = [-0.9, 2.5, 1];
   %B = [1.4, -3, 3.56];      
   
   % Prueba 4. OK
   %E = [-1.9, 2.5, 4];
   %B = [2, 5, 6.56];         
   
   momento_dipolo_mag_ori = 1;
   Dip_E_ori = [0, 0, 0];
   Dip_M_ori = B / norm(B) * momento_dipolo_mag_ori;       
   Dip_M_ori_norm = Dip_M_ori / norm(Dip_M_ori);   
   
   [e1, e2, e3] = ObtenBaseEnh(Dip_M_ori/norm(Dip_M_ori), -v_norm);   
   
   Fmag1_ori = cross(e2, B)
   Fmag2_ori = cross(-e2, B)
   Fmag3_ori = cross(e3, B)
   Fmag4_ori = cross(-e3, B)
   
   Fmag_ori = Fmag1_ori + Fmag2_ori + Fmag3_ori + Fmag4_ori
   
   v1 = Vel_Addition_Law( e2, v_boost);
   v2 = Vel_Addition_Law( -e2, v_boost);
   v3 = Vel_Addition_Law( e3, v_boost);
   v4 = Vel_Addition_Law( -e3, v_boost);
   
   Dopler1 = Sigma * (1 - dot(e2, v_boost));
   Dopler2 = Sigma * (1 - dot(-e2, v_boost));
   Dopler3 = Sigma * (1 - dot(e3, v_boost));
   Dopler4 = Sigma * (1 - dot(-e3, v_boost));      
   
   % Efecto de que la fuerza es sobre partículas en movimiento con
   % diferente dilatació temporal
   F_shit = Fmag1_ori * Dopler1 + Fmag2_ori * Dopler2 + Fmag3_ori * Dopler3 + Fmag4_ori * Dopler4
   
   F_shit_boost = Boost( [0, F_shit], v_boost);
   
   % Cálculo en laboratorio  
   [E_lab, B_lab] = Boost_EM(E, B, v_boost);
   
   Fmag1_lab = cross(v1, B_lab) + E_lab;
   Fmag2_lab = cross(v2, B_lab) + E_lab;
   Fmag3_lab = cross(v3, B_lab) + E_lab;
   Fmag4_lab = cross(v4, B_lab) + E_lab;
   
   F_shit_lab = Fmag1_lab * Dopler1 + Fmag2_lab * Dopler2 + Fmag3_lab * Dopler3 + Fmag4_lab * Dopler4
    
   F_E_ori = 4 * E;
   F4_E_lab = Boost( [0, F_E_ori], v_boost);
   
   F_E_lab = F4_E_lab(2:4);
   
   Dif_tot = F_shit_lab - F_E_lab
   
% Pruebas guarrisimas añadidas a 21/05/2018, tras nivel 9.5
elseif nivel_probaturas == 9.7      
   v_boost = [0.4, 0, 0];      
   vel = norm(v_boost);
   v_norm = v_boost / vel;
   Sigma = 1 / (1 - norm(v_boost)^2)^(1/2);
   
   % Prueba 1. OK
   %E = [1, 0, 0];
   %B = [0, 0, 0];
   
   % Prueba 2. Contracción de la barra es cero, se espera pues Sigma veces
   % más momento. OK
   %E = [0, 1, 0];
   %B = [0, 0, 0];   
   
   % Prueba 3. Contracción de la barra ni cero ni Sigma.   OK
   %E = [3, 1, 2];
   %B = [0, 0, 0];
   
   % Prueba 4. Añadir B y ver que no fastidia nada. OK  
   E = [3, 1, 2];
   B = [-1, 0, 6];   
   
   contract = sqrt((E(1)/Sigma)^2 + E(2)^2 + E(3)^2)/ norm(E);
   
   Dip_E_ori = E / norm(E);
   Dip_M_ori = [0,0,0];
   T_boost = Tensor_boosts(v_boost);   
   T_dip_ori = TensorEM(Dip_E_ori, Dip_M_ori);  
   T_dip_lab_boost = (T_boost * T_dip_ori) * T_boost' ;   
   [Dip_E_lab, Dip_M_lab] = ObtenCampoEMTensor(T_dip_lab_boost);
   
   Val = norm(Dip_E_lab) / (contract*Sigma)     % Debe valer 1
   
end




