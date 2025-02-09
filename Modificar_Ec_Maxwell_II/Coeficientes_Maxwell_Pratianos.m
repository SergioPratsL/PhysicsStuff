function [ dt_E, dt_B ] = Coeficientes_Maxwell_Pratianos( E, B, E_rad_norm, B_rad_norm, a, b )
% Esta función debe calcular la variación del campo electromagnético
% necesaria para que con dicha variación se pueda compensar la variación 
% de energía y momento del campo radiado aisalo de la(s) partículas

% Las variables E y B representan el campo eléctrico TOTAL que hay en un
% punto. 
% La variable a es la ganancia o pérdida de energía del campo radiado de 
% la partícula estudiada (es decir, su campo radiado aislado). Por ganancia
% o pérdida de energía se entiende no la evolución temporal de la energía
% sino dt_u  + div(S) = -a
% Análogamente b es la ganancia o pérdida de momento del campo radiado 
% de la partícula estudiada y se entiende como dt_S - Div(TensMaxwell) = b
% PD: el por qué del -Div() en vez +Div de no lo sé

% Se asume que la corrección de E_vac no puede ser en la dirección de E_rad 
% ni la de B_vac en la de B_rad, se hacen rotaciones puesto que simplifican
% tremendamente las fórmulas, aunque siguen siendo muy difíciles.


% Preparar los coeficientes del sistema de ecuaciones
[e1, e2, e3] = ObtenBaseEnh(E_rad_norm, B_rad_norm);

E_nb = ConvierteNuevaBase(E, e1, e2, e3);
B_nb = ConvierteNuevaBase(B, e1, e2, e3);
b_nb = ConvierteNuevaBase(b, e1, e2, e3);

% Ahora uso fórmulas directamente copiadas de una solución simbólica, por
% ello he de cambiar las varaibles de forma un tanto fea...
E1 = E_nb(1);
E2 = E_nb(2);
E3 = E_nb(3);
B1 = B_nb(1);
B2 = B_nb(2);
B3 = B_nb(3);
b1 = b_nb(1);
b2 = b_nb(2);
b3 = b_nb(3);

E1_vac_nb = 0;
E2_vac_nb = (B1^2*E1*b3 - B1*E2^2*a - B1^2*E2*b3 + E2^2*E3*b2 + E2*E3^2*b3 + B1*B3*E2*b1 + B1*B3*E2*b2 + B1*B3*E3*b3 + B1*E1*E2*a + E1*E2*E3*b1)/(- B1^3*E2 + E1*B1^3 + B1^2*B3*E3 - B1*B3^2*E2 + B1*E2^3 - E1*B1*E2^2 + B1*E2*E3^2 - E1*B3*E2*E3);
E3_vac_nb = (B1^2*E1*b1 - E2^3*b2 + B1^2*E2*b2 + B3^2*E2*b2 + B3^2*E3*b3 - E1*E2^2*b1 - E2^2*E3*b3 + B1*B3*E1*b3 + B1*B3*E3*b1 - B1*E2*E3*a + B3*E1*E2*a)/(- B1^3*E2 + E1*B1^3 + B1^2*B3*E3 - B1*B3^2*E2 + B1*E2^3 - E1*B1*E2^2 + B1*E2*E3^2 - E1*B3*E2*E3);
B1_vac_nb = -(B1^2*B3*b1 + B1^2*B3*b2 + B1*B3^2*b3 + B1^2*E1*a - B1^2*E2*a - B1*E2^2*b3 + B1*E1*E3*b1 + B1*E1*E2*b3 + B1*E2*E3*b2 + B3*E1*E3*b3)/(- B1^3*E2 + E1*B1^3 + B1^2*B3*E3 - B1*B3^2*E2 + B1*E2^3 - E1*B1*E2^2 + B1*E2*E3^2 - E1*B3*E2*E3);
B2_vac_nb = 0;
B3_vac_nb = -(B1^2*E3*a - B1^3*b2 - B1^2*B3*b3 - B1^3*b1 + B1*E2^2*b1 + B1*E2^2*b2 + B1*E3^2*b1 + B3*E3^2*b3 - B1*B3*E2*a + B1*E2*E3*b3 + B3*E2*E3*b2)/(- B1^3*E2 + E1*B1^3 + B1^2*B3*E3 - B1*B3^2*E2 + B1*E2^3 - E1*B1*E2^2 + B1*E2*E3^2 - E1*B3*E2*E3);


E_vac_nb = [E1_vac_nb, E2_vac_nb, E3_vac_nb];
B_vac_nb = [B1_vac_nb, B2_vac_nb, B3_vac_nb];

dt_E = ConvierteBaseOriginal(E_vac_nb, e1, e2, e3);
dt_B = ConvierteBaseOriginal(B_vac_nb, e1, e2, e3);

end

