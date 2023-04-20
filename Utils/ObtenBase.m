function [e1, e2, e3] = ObtenBase(v1)
% Base en la que v1 es la primera componente y las otras dos las obtiene
% gPerpendicular

dims = size(v1);

if dims(1) ~= 1  || dims(2) ~= 3
    disp('v1 mal')
    return
end     


veloc = norm(v1);

if veloc == 0
    e1 = [1, 0, 0];
    e2 = [0, 1, 0];
    e3 = [0, 0, 1];
    return 
end

e1 = v1 / veloc;

[e2, e3] = gPerpendicular(e1);