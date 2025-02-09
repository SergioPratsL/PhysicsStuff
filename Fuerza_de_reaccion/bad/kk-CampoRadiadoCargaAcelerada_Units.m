function [Eo, Bo] = CampoRadiadoCargaAcelerada_Units(R, v, a)
% Al final la única unidad a tener en cuenta es la velocidad de la luz
% Tanto las velocidades como aceleraciones tienen que ser en m/s y m/s^2

dist = norm(R);
%perm_elec = 8,854*10^-12;
%perm_mag = 4 * pi * 10^-7;

%fact_E = 1 / (4*pi*perm_elec);
%fact_B = 4*pi / perm_elec;

c = 299792458;

q = 1;  % La carga puede ser lo que me de la gana (de hecho no la usaré).

v_rel = v / c;

Prod_Vectorial = cross( R, cross((R - dist*v_rel), a));

Eo = (1/c)* 1 / (dist - dot(R,v_rel))^3 * Prod_Vectorial;

% Calculamos B en base a E aplicando el factor 1/c
Bo = (1/c)* cross( R, Eo ) / dist;


end


