function B = ObtenBaseMatrix(v1)
% Base en la que v1 es la primera componente y las otras dos las obtiene
% gPerpendicular

[e1, e2, e3] = ObtenBase(v1);

B = [e1; e2; e3];