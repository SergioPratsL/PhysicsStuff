function [Eo, Bo] = CampoRadCargAcc_FormulaSergio( R, v, dv)
% Campo radiado de de una particula que se mueve a velocidad v
% y experimenta una aceleracion, medida desde el laboratorio de dv/dt
% R indica la distancia

% Este campo no viene dado por la formula de Lienard Wiechert sino que es
% la mia.

% Ya no es neesario porque se arreglo la formula de Lienard-Wiechert y ya
% coincide con el valor obtenido a partir del doble boost!!

dist = norm(R); % tambien es el tiempo (con c=1)
R_dot_v = dot(R, v);

S = dist - R_dot_v;

R_dot_dv = dot(R, dv);


veloc = norm(v);
Sigma = fGamma(veloc);

% Con la base elegida, "v" siempre tiene signo positivo por lo que puede
% usarse "veloc" en su lugar
e3 = v / veloc;
[e1, e2] = gPerpendicular(e3);

% Aceleracion en la nueva base:
dv_new_base = CambioBase(e1, e2, e3, dv);

R_new_base = CambioBase(e1, e2, e3, R);

% Obtener el campo en la base en la que la direcicon de la velocidad es e3
% (en mis formulas siempre es la x)

% Componente de dA/dt
E_e3 = - dist * ( dist * dv_new_base(3) + veloc * R_new_base(1) * dv_new_base(1) + veloc * R_new_base(2) * dv_new_base(2) ) / S^3;
% Componte del gradiente del potencial
E_e3 = E_e3 + ( R_new_base(3) * R_dot_dv / S^3 );


E_e1 = - dv_new_base(1) * (dist^2 - dist * R_dot_v) / S^3;
E_e1 = E_e1 + ( R_new_base(1) * R_dot_dv / S^3) ;


E_e2 = - dv_new_base(2) * (dist^2 - dist * R_dot_v) / S^3;
E_e2 = E_e2 + ( R_new_base(2) * R_dot_dv / S^3) ;


B_e3 = ( R_new_base(2) * dv_new_base(1) - R_new_base(1) * dv_new_base(2) ) / S^2;


B_e1 = - R_new_base(2) * (dist * dv_new_base(3) + veloc * R_new_base(1) * dv_new_base(1) + veloc * R_new_base(2) * dv_new_base(2) ) / S^3;
B_e1 = B_e1 + R_new_base(3) * dv_new_base(2) / S^2;


% En este caso va con el signo cambiado
B_e2 = R_new_base(1) * (dist * dv_new_base(3) + veloc * R_new_base(1) * dv_new_base(1) + veloc * R_new_base(2) * dv_new_base(2) ) / S^3;
B_e2 = B_e2 - R_new_base(3) * dv_new_base(1) / S^2;


Eo = E_e1 * e1 + E_e2 * e2 + E_e3 * e3; 

Bo = B_e1 * e1 + B_e2 * e2 + B_e3 * e3;
