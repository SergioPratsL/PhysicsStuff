
nivel_probaturas = 8;

if nivel_probaturas == 1
% Campos radiados puros
    Ei = [1, 0, 0]
    Bi = [0 , -1, 0]
    %Ei = [1, 1.3, 3]
    %Bi = [0 , 0, 0]     % Prueba campo pi

    v_boost = [0.2, -0.15, 0.19]

    u_ini = norm(Ei)^2 / 2 + norm(Bi)^2 / 2;
    S_ini = cross(Ei, Bi);

    E_S_ini = [u_ini, S_ini]

    E_S_fin_boost = Boost(E_S_ini, v_boost)

    [Ef, Bf] = Boost_EM(Ei, Bi, v_boost);

    u_fin = norm(Ef)^2 / 2 + norm(Bf)^2 / 2;
    S_fin = cross(Ef, Bf);

    E_S_fin = [u_fin, S_fin]

    vel = norm(v_boost)
    Sigma = 1 / sqrt(1-vel^2);

    Ratio1 = E_S_fin_boost(1) / E_S_fin(1)
    Ratio2 = E_S_fin_boost(2) / E_S_fin(2)
    Ratio3 = E_S_fin_boost(3) / E_S_fin(3)
    Ratio4 = E_S_fin_boost(4) / E_S_fin(4)


    % Coeficiente mágico cuadra siempre para un campo radiado
    Coeficiente_magico = 1/(1+v_boost(3)) * 1/Sigma

elseif nivel_probaturas == 2            % Mas alla de la carga que viene de frente fallo.
% Campos inducidos con origen en una partícula.
    Ei = [1, 0, 0]
    Bi = [0 , 0, 0]
    
    %v_boost = [0.3, 0, 0]               % Fue bien con campo B=0
    v_boost = [0.2, -.019, 0.25]
    
    TMaxwell = TensorDeMaxwell(Ei, Bi)
    
    u_ini = norm(Ei)^2 / 2 + norm(Bi)^2 / 2;
    S_ini = TMaxwell(1, :);

    E_S_ini = [u_ini, S_ini]

    E_S_fin_boost = Boost(E_S_ini, v_boost)    
    
    
    [Ef, Bf] = Boost_EM(Ei, Bi, v_boost);

    u_fin = norm(Ef)^2 / 2 + norm(Bf)^2 / 2;
    
    TMaxwell_f = TensorDeMaxwell(Ef, Bf)
    
    v_luz = [1, 1, 0, 0];
    v_luz_f = Boost(v_luz, v_boost);
    vect_new = v_luz(2:4);
    vect_new = vect_new / norm(vect_new);
    
    S_fin = TMaxwell_f(1, :) * vect_new(1) + TMaxwell_f(2, :) * vect_new(2) + TMaxwell_f(3, :) * vect_new(3);
    %S_fin = TMaxwell_f(:, 1) * vect_new(1) + TMaxwell_f(:, 2) * vect_new(2) + TMaxwell_f(:, 3) * vect_new(3);
    %S_fin = S_fin';
    
    %S_fin = vect_new(1) * TMaxwell(1,1) + vect_new(2) * TMaxwell(2,2) + vect_new(3) * TMaxwell(3,3);
    %S_fin = S_fin * vect_new;
    
    
    E_S_fin = [u_fin, S_fin]

    vel = norm(v_boost)
    Sigma = 1 / sqrt(1-vel^2);

    Ratio1 = E_S_fin_boost(1) / E_S_fin(1)
    Ratio2 = E_S_fin_boost(2) / E_S_fin(2)
    Ratio3 = E_S_fin_boost(3) / E_S_fin(3)
    Ratio4 = E_S_fin_boost(4) / E_S_fin(4)


    % Coeficiente mágico cuadra siempre para un campo radiado
    Coeficiente_magico = 1/(1+v_boost(1)) * 1/Sigma

elseif nivel_probaturas == 3                
    Ei = [1, 0, 0]
    Bi = [0 , 0, 0]
    
    v_boost = [-0.34, 0.21, 0.25]    
    Sigma = 1 / sqrt(1 - norm(v_boost)^2);
    
    % Anteriormente llamado Coeficiente_magico :P
    Factor_volumen = 1/(1-dot(v_boost,Ei)/norm(Ei)) * 1 / Sigma;
    % Era precisamente el inverso!!
    Factor_volumen = 1 / Factor_volumen 
    
    % estamos ante un caso simplificado (no B en origen)
    u_ini = norm(Ei)^2 / 2 
    S_ini = [(norm(Ei)^2 / 2), 0, 0];

    E_S_ini = [u_ini, S_ini];

    E_S_fin_boost = Boost(E_S_ini, v_boost);    
    
    % Hardcodeado el campo magnético como (0, 0, 0)!
    [Ef, Bf] = Boost_EM(Ei, [0,0,0], v_boost)
    
    % Dado lo chungo del tema, sólo sacar la energía de momento.
    
    %uf = norm(Ef)^2/2 + norm(Bf)^2/2F
    uf_trampas = norm(Ef)^2/2 - norm(Bf)^2/2;
    
    uf_boost_norm = E_S_fin_boost(1)
    
    uf_trampas_norm = uf_trampas * Factor_volumen

    
% Verificar que [u, S] se transforma como la componente temporal de un
% boost --> Me sale que no :(.
elseif nivel_probaturas == 4

    Ei = [1, 0, 0]
    Bi = [0 , 1, 0]
    
    v_boost = [0.3, 0, 0]    
    Sigma = 1 / sqrt(1 - norm(v_boost)^2);
    
    TMomento_Energia = Tens_Momento_Energia(Ei, Bi)
    
    uf_vect = Boost( TMomento_Energia(1,1:4), v_boost);
    Sxf_vect = Boost( TMomento_Energia(2,1:4), v_boost);
    Syf_vect = Boost( TMomento_Energia(3,1:4), v_boost);
    Szf_vect = Boost( TMomento_Energia(4,1:4), v_boost);
    
    %uf_vect = Boost( TMomento_Energia(1:4, 1)', v_boost);
    %Sxf_vect = Boost( TMomento_Energia(1:4, 2)', v_boost);
    %Syf_vect = Boost( TMomento_Energia(1:4, 3)', v_boost);
    %Szf_vect = Boost( TMomento_Energia(1:4, 4)', v_boost);    
    
    TMomento_Energia_fin_boost = [uf_vect; Sxf_vect; Syf_vect; Szf_vect]
    
    uf = uf_vect(1);
    Sxf = Sxf_vect(1);
    Syf = Syf_vect(1);
    Szf = Szf_vect(1);
    
    u_S_boost = [uf, Sxf, Syf, Szf];
    
    
    [Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
    
    TMomento_Energia_fin_EM = Tens_Momento_Energia(Ef, Bf)
    
    u_S_boost_EM = TMomento_Energia_fin(1, 1:4);
    
    
% Ver si el momento se propaga en la dirección en la cual el tensor de 
% Maxwell apunta precisamente en esa dirección
elseif nivel_probaturas == 5
    Ei = [2, 0, 0]
    Bi = [0 , 3, 0]
    %Ei = [1, -3, 6]
    %Bi = [1.2 , -0-9, 0.5]
    
    v_boost = [0, 0.3, 0]    
    Sigma = 1 / sqrt(1 - norm(v_boost)^2);
    
    Factor_volumen = (1-dot(v_boost,Ei)/norm(Ei)) * Sigma;

    % estamos ante un caso simplificado (no B en origen)
    u_ini = norm(Ei)^2 / 2 + norm(Bi)^2 / 2;
    S_ini = [(norm(Ei)^2 / 2), 0, 0];
    
    TMaxwell_i = TensorDeMaxwell(Ei, Bi)
    p_norm_piloto = [0, 1, 0]
    p_norm_piloto = p_norm_piloto / norm(p_norm_piloto);
    
    p_calculo_nuevo = (TMaxwell_i * p_norm_piloto')'
    p_calculo_nuevo_n = p_calculo_nuevo / norm(p_calculo_nuevo);
    check = dot(p_norm_piloto, p_calculo_nuevo)
    
    
    E_S_ini = [u_ini, S_ini];

    E_S_fin_boost = Boost(E_S_ini, v_boost);
    
    [Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
    
    TMaxwell_f = TensorDeMaxwell(Ef, Bf);
    
    % Resolver el sistema de ecuaciones.
    %A = TMaxwell_f - eye(3);
    %B  = [0, 0, 0]';
    %p = linsolve(A, B)
    
    p_norm_piloto = E_S_fin_boost(2:4);
    %p_norm_piloto = [1, -0.3, 0]
    p_norm_piloto = p_norm_piloto / norm(p_norm_piloto);
    
    p_calculo_nuevo = TMaxwell_f * p_norm_piloto';
    
    %p_calculo_nuevo(2) / p_calculo_nuevo(1)
    
    
% Diagonalizar el tensor de Maxwell para sacar el momento longitudinal
elseif nivel_probaturas == 6        
    %Ei = [2, 0, 0]
    %Bi = [0 , 2, 0]    
    Ei = [1, 0, 0]
    Bi = [0 , 0, 0]    
    
    v_boost = [0, 0.3, 0]    
    Sigma = 1 / sqrt(1 - norm(v_boost)^2);
    
    Factor_volumen = (1-dot(v_boost,Ei)/norm(Ei)) * Sigma;    
    
    TMaxwell_i = TensorDeMaxwell(Ei, Bi);
    
    [V, D] = eig(TMaxwell_i)
    
    p_long_ini = Momento_Longitudinal(TMaxwell_i)
    
    S_ini = cross(Ei, Bi)
    % A falta de nada mejor
    u_ini = norm(Ei)^2 / 2 + norm(Bi)^2 / 2;
    
    p_tot_ini = S_ini + p_long_ini
    
    
    E_S_ini = [u_ini, p_tot_ini];

    E_S_fin_boost = Boost(E_S_ini, v_boost);
    p_fin_boost = E_S_fin_boost(2:4)
    
    [Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
    TMaxwell_f = TensorDeMaxwell(Ef, Bf);    
    
    p_long_fin = Momento_Longitudinal(TMaxwell_f);
    
    S_fin = cross(Ef, Bf)
    
    p_tot_fin = p_long_fin + S_fin
    
    p_tot_fin_norm = p_tot_fin * Factor_volumen

% Pro    
elseif nivel_probaturas == 7        
    Ei = [1, 0, 0]
    Bi = [0 , 0, 0]    
    %Ei = [1, 0.5, -1.9];
    %Bi = [1.4 , 2, 0];    
    
    v_boost = [0, 0.4, 0]      % Cuadra con Ei y con campo doble
    %v_boost = [-0.4, 0.3, 0.1]  % Cuadra con Ei  y con campo doblo
    Sigma = 1 / sqrt(1 - norm(v_boost)^2);    
    
    j_i = [1,0,0,0];
    j_f = [1, v_boost] * Sigma;
    
    TMaxwell_i = Tens_Momento_Energia(Ei, Bi);
    
    pu_i = j_i * TMaxwell_i
    
    pf_boost = Boost(pu_i, v_boost)
    
    [Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
    TMaxwell_f = Tens_Momento_Energia(Ef, Bf);
    
    pf = j_f * TMaxwell_f 
    

% Ver si lo que he hecho para momentos dipolares sirve aquí (allí de
% momento no va demasiado bien).
elseif nivel_probaturas == 8
    
    Ei = [1, 10.5, -1.9];
    Bi = [1.4 , 2, 5];    
    
    v_boost = [-0.39, 0.3, 0.25];      % Cuadra con Ei y con campo doble
    Sigma = 1 / sqrt(1 - norm(v_boost)^2);        
    
    T_boost = Tensor_boosts(v_boost);
    
    T_Mom_E_i = Tens_Momento_Energia(Ei, Bi);
    
    T_Mom_E_boost = (T_boost * ((T_boost * T_Mom_E_i)'))'
    
    [Ef, Bf] = Boost_EM(Ei, Bi, v_boost);
    
    T_Mom_E_direct = Tens_Momento_Energia(Ef, Bf)
    
    
    
end
    
