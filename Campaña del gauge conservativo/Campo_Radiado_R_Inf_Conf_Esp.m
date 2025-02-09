Campo_Radiado_R_Inf_Conf_Esp

% Este script busca comprobar una hipotesis que formule para analizar como
% se puede compensar la energia del campo radiado.
% La prueba considera que dos particulas C y D estan separadas una distancia L en
% el eje x. Estamos en el SRI en el que temporalmente esta C y D se mueve
% con velocidad vx. La clave es que la posicion y el movimiento ambos son
% sobre el eje x.

% Hay que evaluar el campo en diferentes angulos. Al final solo importa el
% coseno entre el vector R y el eje x ya que a diferentes cosenos, le
% corresponde un momento diferente de la radiacion de D.

% Masa y carga normalizadas a 1, c obviamente tambien :P
% A pesar de ser infinito, a la distancia le quitare el factor 1/R!!

% Siempre eje x!
D = 10;     % posicion actual en eje x de D respecto a C
v = -0.4;

% Angulo siempre sobre el eje XY, es indiferente que sea Y o Z.
Rho = pi/4;     % Al ser R infinito, el angulo sera el mismo desde C que desde D

%R lo normalizo a 100, pero daria igual si lo hiciera a 1.
R = 100 * [cos(Rho), sin(Rho), 0];

D_ret = D / (1+v);

R_aux = [D_ret, 0, 0];
v_d = [v, 0, 0];

% Fuerza que sufre C "ahora"
[Ec,H] = CampoInducidoLienardWiechert( -R_aux, v_d );

Fc = E;
ac = Fc;

dv_c = [ac, 0, 0]; 
v_c = [0, 0, 0];

[E_rad_c, H_rad_c] = CampoRadiadoCargaAcelerada( R, v_c, dv_c)


% Momomento relevante de D.
Desfase = D * cos(Rho);     % Desfase en tiempo
D_desf = D + Desfase * v;

% Fuerza que sufre D en ese momento






