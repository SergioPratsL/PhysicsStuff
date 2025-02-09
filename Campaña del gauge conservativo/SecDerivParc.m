function [ddV, ddAx, ddAy, ddAz] = SecDerivParc(Rx, Ry, v)

% Funcion que te saca la segunda derivada para minimizar el riesgo de
% erratas en los calculos

% Solo obtendre los valores que vaya a usar y el resto los dejo a infinito
% para que si hay un valor que deberia tener pero no tengo lo sepa!

R = sqrt( Rx^2 + Ry^2);

ddV = [inf, inf, inf; inf, inf, inf; inf, inf, inf];

ddAx = [inf, inf, inf; inf, inf, inf; inf, inf, inf];

ddAz = [inf, inf, inf; inf, inf, inf; inf, inf, inf];

ddAy = [inf, inf, inf; inf, inf, inf; inf, inf, inf];


ddV(1,1) = (2*Rx^2 - Ry^2) / R^5;
ddV(2,2) = (2*Ry^2 - Rx^2) / R^5;

ddAx(1,1) = (2*Rx^2 - Ry^2) / (v*R^5);
ddAx(2,2) = (2*Ry^2 - Rx^2) / (v*R^5);

ddAy(1,1) = 3 * Ry * Rx / (R^5*v);
ddAy(2,2) = (-Rx/v) * (2/(Ry^3*R) + 1/(Ry*R^3) + 3*Ry/R^5) - 2/(abs(v)*Ry^3);



ddV(1,2) = (3 * Rx *  Ry) / R^5;
ddV(2,1) = ddV(1,2);

ddAx(1,2) = (3 * Rx*  Ry) / (v*R^5);
ddAx(2,1) = ddAx(1,2);

ddAy(1,2) = (2*Ry^2 - Rx^2) / (v*R^5);
ddAy(2,1) = ddAy(1,2);