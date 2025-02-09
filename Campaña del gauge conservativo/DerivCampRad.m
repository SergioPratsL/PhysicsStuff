function [dE, dB] = DerivCampRad( Rx, Ry, a, dR )

% Esta funcion da las derivadas del campo para una variacion de posicion
% en un contexto muy especifico...

% · Las dos particulas (emisora y "punto" receptor estan en el plano XY
% · La particula emisora no se mueve

% a se emplea para sacar el campo radiado y dR para determinar que
% derivadas queremos

% Rx y Ry es la distancia pasada del emisor, no la actual!!

R = [Rx, Ry, 0];

dist = norm(R);

R_norm = R / dist;

dR_norm = dR / dist;

Prod_esc = dot(R_norm, dR_norm);


R_ort = cross(R_norm, a);

R_rad = cross(R_norm, R_ort);

% La derivada total sera:
% dE = E_num * dE_Denom + dE_num * E_Denom

% El denominador es 1/R^3 por lo que crecera como 1/R^3*(1-3dr·(R/|R|))
% Pasare a usar 1/R y no 1/R^3...
dE_Denom = - dot(R_norm, dR_norm) / dist^2; 
dB_Denom = dE_Denom;

E_Denom = 1 / dist;
B_Denom = E_Denom;

E_num = R_rad;
B_num = cross( R_norm, R_rad );

% 10.01.2016
v_0 = [0, 0, 0];
[E, B] = CampoRadiadoCargaAcelerada(R, v_0, a);

% La diferencia en el numerador es el resultado de hacer
% (R + dR) x ( (R + dR ) x a)

% Hay que tener en cuenta que dR es infinitesimal por lo que hay que
% descartar los termnos con dos dRs!!

% Por tanto tenemos:  (R + dR) x ( (R + dR) x a )

% Aplicando las identidades A x (B x C) = (A·C)B - (A·B)C da


% dE_num = (R_norm·a)dR_norm + (dR_norm·a)R_norm - 2*(dR_norm·R_norm)a
% dE_num = dot(R_norm, a) * dR_norm + dot(dR_norm, a) * R_norm - 2*dot(R_norm, dR_norm)* a;  

% No exactamente, en cada caso al hacer (R_norm + dR_norm) hay que
% renormalizar por si dR · R <> 0...
% En vez de (R_norm + dv_norm) debe haber (R_norm + dv_norm)*( 1 - R_norm·dv_norm)

% 03.01.2016: ¿No sera mas bien 
% (R_norm + dv_norm)*( 1 - 2*R_norm·dv_norm)??? ... No...


% Tras algo de algebra, lo que queda es:
% (1- 2*R_norm·dR_norm)* [(R_norm·a + dR_norm·a)(R_norm + dR_norm) - (R^2+ 2*R·dr)a]

% Los terminos nuevos de orden 1 son:
% - 2*R_norm·dR_norm* [(R_norm·a)R_norm - a] + (R_norm·a)dR_norm + (dR_norm·a)R_norm - 2*(dR_norm·R_norm)a

% Notese que hay dos terminos que se eliminan :)

dE_num = dot(R_norm, a) * dR_norm + dot(dR_norm, a) * R_norm - 2*dot(R_norm, dR_norm)* dot(R_norm,a) * a;  

% Y para el campo magnetico hay que sacar (R_norm + dR_norm) (1- R_norm·dR_norm) x (E_num + dE_num)

%dB_num = cross(R_norm, dE_num) + cross(dR_norm, E_num);
dB_num = cross(R_norm, dE_num) + cross(dR_norm, E_num) - dot(R_norm, dR_norm)*cross(R_norm, E_num);

% Ni un paso atras, ni un paso en falso!!!
%dB_num = cross(R_norm, dE_num); % Blasfemia

% Juntamos los terminos
% Cambio 10.01.2016. Uso el old para comparar los dos.
dE_old = dE_Denom * E_num + dE_num * E_Denom;

dE = dot(dR_norm,R_norm)*a + dot(R_norm, a)*dR_norm - 3*dot(dR_norm, R_norm)*dot(R_norm, a)*R_norm + dot(dR_norm,a)*R_norm;
dE = dE / dist;

% Cambio 10.01.2016
dB_old = dB_Denom * B_num + dB_num * B_Denom;

dB = -dot(dR_norm,R_norm)*B + cross(dR_norm, E) + cross(R_norm, dE);

check_E = dot(dE,E);
check_B = dot(dB,B);

check = check_E / check_B;



