
% Rx_rep y Ry_rep es el vector retardado visto desde el receptor (plano ZY)
% v es la velocidad con la que el receptor ve al emisor...

% Caso 1. Success!
% v = -0.4;
% Rx_rep = v * Sigma;
% Ry_rep = 1;

% Caso 2. Success!!!
% Rx_rep = 2;
% Ry_rep = 1;
% v = 0.3;


% Caso 2. Success!!!
Rx_rep = 1.5;
Ry_rep = 4;
v = 0.42;

Sigma = 1 / sqrt( 1 - v^2);

% Obtener el vector que ve el emisor
t_rep = sqrt( Rx_rep^2 + Ry_rep^2);
V_rep = [t_rep, Rx_rep, Ry_rep, 0];
v_vect = [v, 0, 0];


R_emis = Boost( V_rep, v_vect);

% -v porque toda formula la hice desde el punto de vista del emisor
dx = DerivadasPotencialYinn( R_emis(2), R_emis(3), -v, 'x');
dy = DerivadasPotencialYinn( R_emis(2), R_emis(3), -v, 'y');
dz = DerivadasPotencialYinn( R_emis(2), R_emis(3), -v, 'z');
dt = DerivadasPotencialYinn( R_emis(2), R_emis(3), -v, 't');

% dz_rep = dz;
% dy_rep = dy;
dy_rep = Boost( dy, -v_vect );
dz_rep = Boost( dz, -v_vect );

% Vuelvo al sistema del receptor los 4 vectores recibidos
dx_aux = Boost( dx, -v_vect );
dt_aux = Boost( dt, -v_vect );

% Transformo las coordenadas de la derivada
% ...vigliar el signo.
dx_rep = Sigma * ( dx_aux - v * dt_aux );
dt_rep = Sigma * ( dt_aux - v * dx_aux );


Ex = - dx_rep(1) - dt_rep(2);
Ey = - dy_rep(1) - dt_rep(3);
Ez = - dz_rep(1) - dt_rep(4);

Bx = dy_rep(4) - dz_rep(3);
By = dz_rep(2) - dx_rep(4);
Bz = dx_rep(3) - dy_rep(2);

E_yinn = [Ex, Ey, Ez]
B_yinn = [Bx, By, Bz]

% Comparar con el campo de L-W...

R = [Rx_rep, Ry_rep, 0];
[E_LW, B_LW] = CampoInducido_sin_unidades(R, v_vect)