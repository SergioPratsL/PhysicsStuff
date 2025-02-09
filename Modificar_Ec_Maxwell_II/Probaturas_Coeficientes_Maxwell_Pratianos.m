

% Las primeras probaturas servir�n para tener una idea de los coeficientes
% obtenidos alrededor de una carga esf�rica quieta donde el campo 
% dominante es el de la carga.

%q = 0.001;  % Valor de la carga, lo pongo tan peque�o para que el campo acelerado no barra al campo externo
q = 1;
c = 1;      % Que no quede duda que la velocidad de la luz es 1
r = 0.001;  % esta cutrez representa la "esfera peque�ita"
m = 0.27;  % Para no ser tan aburrido pongo esto en vez de 1


% Caso 1: carga "casi puntual" radiando, nos colocamos en +X.
% los valores iniciales ser�n un poco guarros:

E_ext = [0, 4, 0]
B_ext = [0, 0, 1]

v = [0, 0, 0];

acel_vect = FLorentz(v, E_ext, B_ext) * q / m;


vector_direccion_unitario = [1, 0, 0];
pos_evaluada = r * vector_direccion_unitario;

[E_rad, B_rad] = CampoRadiadoCargaAcelerada( pos_evaluada, v, acel_vect)

E_rad_norm = E_rad / norm(E_rad); 
B_rad_norm = B_rad / norm(B_rad);

acel = norm(acel_vect);       

% Representa la densidad de energ�a radiada en r, con un mont�n de
% constantes olvidadas.
a = (acel*q/r)^2;
% Representa el momento radiado, si estamos mirando en la direcci�n +X
b = a * vector_direccion_unitario;

% Campo radiado por la carga
E_int = [1, 0, 0] / r^2
% La carga est� quieta y si hubiere campo magn�tico externo, palidece en
% comparaci�n con el campo el�ctrico interno.
B_int = [0, 0, 0]          

E = E_ext + E_rad + E_int
B = B_ext + B_rad + B_int

%[dt_E, dt_B] = Coeficientes_Maxwell_Pratianos_bad1( E, B, a, b )
[dt_E, dt_B] = Coeficientes_Maxwell_Pratianos( E, B, E_rad_norm, B_rad_norm, a, b )

