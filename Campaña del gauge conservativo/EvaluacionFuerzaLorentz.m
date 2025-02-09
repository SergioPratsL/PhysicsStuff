

% Prueba 1 (coinciden)
%v = [0.3, 0, 0];
% nuevo SRI
%u = [0, 0.3, 0];
%E = [0, 1, 0];
%B = [0, 0.1, 0.2];

% Prueba 2 (coinciden)
v = [0, 0, 0];
% nuevo SRI
u = [0, 0.3, 0];

E = [0, 1, 1];
B = [0.1, 0.1, 0.2];


% Prueba 3
v = [0.2, 0, -01];
% nuevo SRI
u = [0.1, -0.3, 0.18];

E = [0, 1, 1];
B = [0.2, 0.1, -0.15];


F =  [0, 0, 0]; % inicializar

F = fLorentz( v, E, B);

Pot = dot(F, v);

Momento_origen = [Pot, F(1), F(2), F(3)]

Momento = Boost( Momento_origen, u)
 
% Transformar el intervalo de tiempo y sacar la nueva velocidad
vel_4 = [1, v(1), v(2), v(3)];
vel_bost_4 = Boost(vel_4, u); 
t_bis = vel_bost_4(1);

v_bis = [vel_bost_4(2), vel_bost_4(3), vel_bost_4(4)] / t_bis;

% Transformar el campo de SRI
[E_bis, B_bis] = Boost_EM(E, B, u);

F_bis = fLorentz( v_bis, E_bis, B_bis );

Pot_bis = dot(F_bis, v_bis);

Momento_bis = [Pot_bis, F_bis(1), F_bis(2), F_bis(3)] * t_bis
