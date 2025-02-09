

% DECLARACION DE CONSTANTES en el Sistema de Referencia Internacional:

% Radio de Bohr:   5.2917721067(12)×10?11 m
a0 = 5.2918 * 10^-11;
% Carga del electron:  1.602 176 565(35)×10?19  Coulombs
e = 1.6022 * 10^-19;

h = 6.6261 * 10^-34;
h_bar = h / (2*pi)

me = 9.11 * 10^-31 


Permitividad_vacio = 8.85419 *10^-12;       % Vergonzoso, me lo habia dejado


factor1 = h_bar^2 / (2* me * e * a0^2)

% (1/e) para expresar en eV
factor2 = e^2 / (4*pi*Permitividad_vacio*a0) * (1/e) 


syms r;
syms a;

Ond1 = exp(-r/a);
d_Ond1 = diff(Ond1, r);

r2_d_Ond1 = r^2 * d_Ond1;

d_r2_d_Ond1 = diff(r2_d_Ond1, r);

% Ahora vamos para la integral, a ver, la última parte de la derivación
% esférica implica un 1/r^2 pero integrar implica un 4pi*r^2 así que al
% final queda el 4pi.

pre_int_Ond1 = 4 * pi * d_r2_d_Ond1 * Ond1;

int_Ond1 = int(pre_int_Ond1)
% int_Ond1 = (pi*exp(-(2*r)/a)*(a^2 + 2*a*r - 2*r^2))/a

val_1 = - subs(int_Ond1, r, 0)
% pi * a


Ond2 = (2 - r/a) * exp(-r/(2*a));
d_Ond2 = diff(Ond2, r);

r2_d_Ond2 = r^2 * d_Ond2;

d_r2_d_Ond2 = diff(r2_d_Ond2, r);

pre_int_Ond2 = 4 * pi * d_r2_d_Ond2 * Ond2;

int_Ond2 = int(pre_int_Ond2)
% int_Ond2 = (pi*exp(-r/a)*(8*a^4 + 8*a^3*r - 12*a^2*r^2 + 8*a*r^3 - r^4))/a^3

val_2 = - subs(int_Ond2, r, 0)


Ond1_quad = Ond1 * Ond1;

pre_int_pot_Ond1 = 4 * pi * r * Ond1_quad;

int_Ond1_pot = int(pre_int_pot_Ond1)

val_1_pot = - subs(int_Ond1_pot, r, 0)



Ond2_quad = Ond2 * Ond2;

pre_int_pot_Ond2 = 4 * pi * r * Ond2_quad;

int_Ond2_pot = int(pre_int_pot_Ond2)
% int_Ond2_pot = -(4*pi*exp(-r/a)*(2*a^3 + 2*a^2*r - a*r^2 + r^3))/a

val_2_pot = - subs(int_Ond2_pot, r, 0)

% El resultado es que  la corrección no es necesaria y nada de lo que he
% hecho vale una mierda :@!!   ... pas si vite!!!



pre_int_Ond1_quad = 4 * pi * r^2 * Ond1_quad;

%q_int_Ond1 = int(pre_int_Ond1_quad)
%test2 = diff(pre_int_Ond1_quad, r)     %% mal!

q_int_Ond1 = int(pre_int_Ond1_quad, 0, r)
q1_en_0 = subs(q_int_Ond1, r, 0)    %% ok, da 0


E1 = q_int_Ond1 / r^2;
%xxx = eval(subs(E1, {r, a}, {a0/1000, a0}));

Ep_E1 = (1/r^2) * E1;

pre_int_Ep_E1 = 4 * pi * r^2 * Ep_E1;

int_Ep_E1 = int(pre_int_Ep_E1)

%val_energ_Ep_E1 = - subs(int_Ep_E1, r, 0)
val_energ_Ep_E1 = - eval(subs(int_Ep_E1, {r, a}, {a0/1000, a0}))

energ_EpE1 = e^2/(4*pi*Permitividad_vacio)^2 * (Permitividad_vacio/2) * 1/(pi*a0^3) * val_energ_Ep_E1 * (1/e)


pre_int_Ond2_quad = 4 * pi * r^2 * Ond2_quad;

q_int_Ond2 = int(pre_int_Ond2_quad, 0, r)
q2_en_0 = subs(q_int_Ond2, r, 0)    %% ok, da 0


E2 = q_int_Ond2 / r^2;

Ep_E2 = (1/r^2) * E2;

pre_int_Ep_E2 = 4 * pi * r^2 * Ep_E2;

int_Ep_E2 = int(pre_int_Ep_E2)

val_energ_Ep_E2 = - eval(subs(int_Ep_E2, {r, a}, {a0/1000, a0}))

energ_EpE2 = e^2/(4*pi*Permitividad_vacio)^2 * (Permitividad_vacio/2) * 1/(32*pi*a0^3) * val_energ_Ep_E2 * (1/e)



% Los valores que habrían si no hubiera corrección de la energía
E1_quad = E1 * E1;
pre_int_E1_quad = 4 * pi * r^2 * E1_quad; 
int_E1_quad = int(pre_int_E1_quad)
val_energ_E1_quad = - eval(subs(int_E1_quad, {r, a}, {a0/1000, a0}))

energ_EpE1 = e^2/(4*pi*Permitividad_vacio)^2 * (Permitividad_vacio/2) * 1/(pi*a0^3)^2 * val_energ_E1_quad * (1/e)


E2_quad = E2 * E2;
pre_int_E2_quad = 4 * pi * r^2 * E2_quad; 
int_E2_quad = int(pre_int_E2_quad)
val_energ_E2_quad = - eval(subs(int_E2_quad, {r, a}, {a0/1000, a0}))

energ_EpE2 = e^2/(4*pi*Permitividad_vacio)^2 * (Permitividad_vacio/2) * 1/(32*pi*a0^3)^2 * val_energ_E2_quad * (1/e)