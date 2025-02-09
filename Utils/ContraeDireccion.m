
function R2 = ContraeDireccion(R1, v)
% Contrae Sigma la dirección de V, respeta el resto de direcciones.

veloc = norm(v);

if veloc >= 1
    disp('mal')
    return
end 
    
if veloc == 0
    R2 = R1;
    return 
end 

Sigma = fGamma(veloc);

e_v = v / veloc;

prod_esc = dot(R1, e_v);

R2 = R1 - (1-1/Sigma)*prod_esc * e_v;


