function [e1, e2, e3] = ObtenBaseEnh( v1, v2 )
% Base en la que v1 es la primera componente, v2 es la segunda componente
% "favorita" y v_3 la restante.
% Se han de pasar vectores normalizados

if norm(cross(v1,v2)) == 0
    [e1, e2, e3] = ObtenBase(v1);
    return
end

e1  = v1;
e2 = v2 - dot(v1, v2) * v1;
e2 = e2 / norm(e2);
e3 = cross(e1, e2);


end

