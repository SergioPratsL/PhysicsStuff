
function Vout = CambioBase(e1, e2, e3, Vin)
% e1, e2 y e3 deben ser ortonormales y de dimension 3

Vout(1) = dot(Vin, e1);
Vout(2) = dot(Vin, e2);
Vout(3) = dot(Vin, e3);
