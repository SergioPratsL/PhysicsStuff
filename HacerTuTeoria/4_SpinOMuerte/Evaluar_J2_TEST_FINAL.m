% Script para evaluar el valor de J^2 en diferentes direcciones para un
% estado con j=3/2, jz=1/2, si no me equivoco, este estado debe una CL
% de los estados (2,1,0,+) y (2,1,1,-) ambos con amplitud 1/raiz(2)

% Este Matlab sacará el resultado del término L·S sobre esta combinación 
% de estados para ver si es proporcional al estado origen y por tanto
% cumple estrictamente la condición de autovector...

S_mas = [1,0];
S_menos = [0,1];

% Probaré valores arbitrarios de Theta y Fi
Theta = pi/2;
Fi = 0;
 
%Theta = 0;
%Fi = pi/4;


%Theta = pi/2;
%Fi = pi/4;


%Theta = pi/4;
%Fi = pi/3;


%Theta = 0.7054 * pi;
%Fi = 0.295 * pi;


%Theta = 1.45 * pi;
%Fi = -0.398 * pi;


%Theta = 0.29 * pi;
%Fi = -0.318 * pi;

raiz_2 = sqrt(2);
 
% Valores buenos (no cambiar!)
cof_1 = sqrt(2);
%cof_1 = -1/sqrt(2);
cof_2 = 1;          % Siempre ha de valer 1 con mi enfoque.

% https://en.wikipedia.org/wiki/Spherical_coordinate_system#Integration_and_differentiation_in_spherical_coordinates
dir_fi = [-sin(Fi), cos(Fi), 0]; 
dir_theta = [cos(Fi)*cos(Theta), sin(Fi)*cos(Theta), -sin(Theta)]; 

Px = MatrizPauli(1);
Py = MatrizPauli(2);
Pz = MatrizPauli(3);

% Matrices en las direcciones ortogonales de Theta y de Fi
P_theta = Px * dir_theta(1) + Py * dir_theta(2) + Pz * dir_theta(3);
P_fi = Px * dir_fi(1) + Py * dir_fi(2) + Pz * dir_fi(3); 

% Calentanmiento 1: Sacar Lz en el estado |1,1>, sin spin
% Esto ha ido bien a la primera... La ratio siempre es 1
%  f = exp(1i*Fi)*sin(Theta)

%  % Debido a que dir_fi no tiene dirección en Z, podemos ignorarla.
%  Lz = exp(1i*Fi)*(-1i)^2*dir_theta(3)
%  
%  ratio_tranining_1 = Lz / f
 
% Calentamiento 2: Sacar Lz sobre el estado |1,1,->
% También fue bien a la primera.
% f = exp(1i*Fi)*sin(Theta)*S_menos'
% 
% Lz = exp(1i*Fi)*(-1i)^2*dir_theta(3)*S_menos'+1/2*Pz*f
% 
% ratio_tranining_2 = Lz(2) / f(2)


% El primer raiz de 2 es porque las ondas del átomo H no tienen igual
% amplitud, es más "grande" la (2,1,0)
% Hay que añadir otro raiz(2) puesto que la proporción no era igual para
% los dos estados sino que ms=1/2 era raiz(2/3) y ms=-1/2 era raiz(1/3)
% Aunque parece absurdo decir raiz_2^2, es para indicar sus dos orígenes

% ¡¡OJO!! ESTAS FUNCIONES HAN SIDO DEFINIDAS CON Theta siendo 0 en el polo
% sur y pi en el norte, exactamente lo contrario que lo que se usa en la
% wikikedia, de allí que voy a ponerr unos bonitos signos menos, jeje
% NOTA: me estaba volviendo loco con el infierno polar... esto me recuerda
% a "en las montañas de la locura", que también ocurre en los polos
% Lo menos habré perdido 1D6 de cordura por intentar sacar esto
F1 = cof_1*raiz_2*cos(-Theta)*S_mas';
F2 = exp(1i*Fi)*sin(-Theta)*S_menos';

f = F1 + F2

%f = cof_1*raiz_2*cos(Theta)*S_mas + cof_2*exp(1i*Fi)*sin(Theta)*S_menos

%Descomponerlo al máximo:
T1 = cof_1*raiz_2*sin(-Theta)*S_mas';
T2 = -exp(1i*Fi)*cos(-Theta)*S_menos';
T3 = 1i*exp(1i*Fi)*S_menos';

L_wave_vect_x = Px*((T1+T2)*dir_fi(1)+T3*dir_theta(1));
L_wave_vect_y = Py*((T1+T2)*dir_fi(2)+T3*dir_theta(2));
% %L_wave_vect_z = Pz*((T1+T2)*dir_fi(3)+T3*dir_theta(3));
L_wave_vect_z = Pz*(T3*dir_theta(3));       % Equivalente

L_dot_S = L_wave_vect_x + L_wave_vect_y + L_wave_vect_z;

% factor -i que relaciona el momento y el gradiente
L_dot_S = - 1i * L_dot_S

%L_dot_S_original = cof_1*raiz_2*(-sin(Theta))*P_fi*S_mas' + cof_2*exp(1i*Fi)*(cos(Theta)*P_fi - 1i*P_theta)*S_menos';

ratio1 = L_dot_S(1) /f(1)
ratio2 = L_dot_S(2) /f(2)
ratio_total = norm(L_dot_S) / norm(f);




% KKs
X =  1 / sqrt(1+(1+1/raiz_2)^2);
Y = 1 / sqrt(1+((1+sqrt(5))/2)^2);

B = -1;
C = -1;
A = 1;
Answer1= ((-B)+((B^2-4*A*C))^0.5)/(2*A);
Answer2= ((-B)-((B^2-4*A*C))^0.5)/(2*A);