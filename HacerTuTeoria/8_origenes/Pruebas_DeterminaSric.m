% Probar el algoritmo que creó Gemini.


%E = [1, 0, 0];
%B = [0, 0.1, 0];
% v_sric = 0         0    0.1000

%E = [1, 0, 0];
%B = [-0.24, 0.1, 0.4];
% v_sric = 0   -0.3744    0.0936

% En este caso B gana a E
%E = [-0.2, 0.3, 0];
%B = [0.9, 1.1, 0.4];
% v_sric = 0.0548    0.0365   -0.2237

% E = [1, 0, 0];
% B = [0, 0, 0];
% v_sric =  0     0     0


E = [1, 0, 0];
B = [0, 1, 0];
% v_sric =  0     0     1



S = cross(E, B)
v_sric = DeterminaSric(E, B)
[E_prim, B_prim] = Boost_EM(E, B, v_sric);
S_prim = cross(E_prim, B_prim)


