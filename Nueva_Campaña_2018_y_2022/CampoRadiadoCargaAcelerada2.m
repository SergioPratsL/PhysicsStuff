function [Eo, Bo] = CampoRadiadoCargaAcelerada2(R, v, a)
% Normalizado como siempre, creo esta funci�n porque no me f�o de la
% anterior

dist = norm(R);

Prod_Vectorial = cross( R, cross((R - dist*v), a));

Eo = 1 / (dist - dot(R,v))^3 * Prod_Vectorial;

Bo = cross( R, Eo ) / dist;


end

