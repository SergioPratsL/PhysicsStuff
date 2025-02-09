clear

% Estas ecuaciones asumen que la direccion e1 es la propagacion del campo
% electrico radiado, e2 es la del campo magnetico radiado y e3 la de la
% propagacion y asumo la condicion que la variacion del campo causada por
% la correccion (o por el vacio) es nula en la direccion del campo radiado
% por lo que E1_vac = 0 y B2_vac = 0

E1 = sym('E1');
E2 = sym('E2');
E3 = sym('E3');

B1 = sym('B1');
B2 = sym('B1');
B3 = sym('B3');

E1_vac = sym('E1_vac');
E2_vac = sym('E2_vac');
E3_vac = sym('E3_vac');
B1_vac = sym('B1_vac');
B2_vac = sym('B2_vac');
B3_vac = sym('B3_vac');

E1_rad = sym('E1_rad');
E2_rad = sym('E2_rad');
E3_rad = sym('E3_rad');

B1_rad = sym('B1_rad');
B2_rad = sym('B2_rad');
B3_rad = sym('B3_rad');


a = sym('a');
b1 = sym('b1');
b2 = sym('b2');
b3 = sym('b3');

% 03.05: el factor 2*a fue un fallo garrafal!
%eqn1 = E1_vac*E1 + E2_vac*E2 + E3_vac*E3 + B1_vac*B1 + B2_vac*B2 + B3_vac*B3 == -2*a;
eqn1 = E1_vac*E1 + E2_vac*E2 + E3_vac*E3 + B1_vac*B1 + B2_vac*B2 + B3_vac*B3 == -a;
eqn2 = E2_vac*B3 - E3_vac*B2 - B2_vac*E3 + B3_vac*E2 == -b1;
eqn3 = -E1_vac*B3 + E3_vac*B1 + B1_vac*E3 - B3_vac*E1 == -b2;
eqn4 = E1_vac*B2 - E2_vac*B1 - B1_vac*E2 + B2_vac*E1 == -b3;
%eqn5 = E1_vac*E1_rad + E2_vac*E2_rad + E3_vac*E3_rad == 0;
%eqn6 = B1_vac*B1_rad + B2_vac*B2_rad + B3_vac*B3_rad == 0;
eqn5 = E1_vac == 0;
eqn6 = B2_vac == 0;

[Output1, Output2] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6], [E1_vac, E2_vac, E3_vac, B1_vac, B2_vac, B3_vac])

Result = linsolve(Output1, Output2)



