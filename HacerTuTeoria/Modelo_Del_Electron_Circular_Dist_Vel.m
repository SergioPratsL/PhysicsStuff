% Este modelo sólo tiene dos elementos pero la gracia es que uno va a
% diferente velocidad

% Resolver 2v^3-v^2+1=0
% Resolver 2v^3-v^2+a=0

%a = 1;  % m1 =   0.6983,  q2 = 0 :(
%a = -1; % requiere v > 1!!

% a = 1/2; % v = -0.8324  m1 = 0.5907  q1 = -0.7271   L = 1/4

%a = 1.5; % v = -0.4554 m1 = 0.8282 L = 3/4. Este es el que más me gusta.

%a = 2;  % Cargas infinitas y M NaN :(.

%a = 1.99; % Todas las raices de las velocidades son positivas...

%a = 1.7; % Aún una sola raiz real para v (-0.3531) y masas correctas.

%a = 1.78815; % Por encima de esto 3 soluciones reales para la velocidad...
%cosa mala. Aquí m1 = 0.9184   v = -0.2981

%a = 1.75; % v=-0.323   m1=0.9055 

%a = 1.7882; % Con 1.7883 ya aparecen 3 soluciones.

%a = 5.105;  % Cuadra fuerza vs energía, pero requiere ir más rápido que la luz. v = 2.7568

% Todo a negativo requiere superar la velocidad de la luz.
%a = -0.1 %a = -1; %a = -3;

a = 1.01;

p = [2,-a,0,2-a];
sols = roots(p);

if norm(imag(sols(1))) < 10^-8
    v = sols(1)
end

if norm(imag(sols(2))) < 10^-8
    v = sols(2)
end

if norm(imag(sols(3))) < 10^-8
    v = sols(3)
end

if abs(v) > 1
    disp('Supera la velocidad de la luz!')
    return
end

m1 = 1 / (1+v^2)

m2 = 1 - m1

% Parad las trampas al solitario, para la carga es imposible porque 
% q1 + q2 = -1
% q1 + vq2 = -1
% No se puede cambiar

L = m1 + m2 * v

q1 = - (a-v)/(1-v)
q2 = - (1-a)/(1-v)

M = q1 + v*q2


me = 9.1093 * 10^-31;
c = 299792458;
h = 6.62607 * 10^-34;
h_bar = h /(2*pi);
q = 1.6022 / 10^(19);
perm_elec = 8.854 *  10^-12;

E_elec = me * c^2;

% ¡Recordar que m1, m2 y v estan todos normalizados!
r = h_bar/2 * 1/(c*me*(m1+m2*v))

% Longitud de onda de compton.
r_compton = h/(c*me);

r_div_r_compton = r / r_compton


% Otras cosas
r_bohr =  5.291 * 10^-11;
r_bohr_div_r_compton = r_bohr / r_compton;

% Freceuencia a la que gira la componente 1
w1 = c / r;

% Esto es la fuerza a una distancia obtenida por otro medio!
F_ideal = -1/(4*pi*perm_elec)*q^2*q1*q2/r^2;

acel_ideal = F / (me*m1)
acel_necesaria = c^2 / r



% FALLÓ.

 




