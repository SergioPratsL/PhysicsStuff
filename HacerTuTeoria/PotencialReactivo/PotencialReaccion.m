function PR = PotencialReaccion(v, a_prop)
% ESTA ES LA FUNCIÓN QUE SE QUIERE PROBAR QUE ES VÁLIDA EN ESTA CAMPAÑA

% v es la velocidad a la que vemos que se mueve el objeto y 
% a_prop es la aceleración medida en el sistema propio del SRI.

%PR_prop = (2/3) * [-dot(a_prop,a_prop), a_prop];
PR_prop = (2/3) * [0, a_prop];

PR = Boost(PR_prop, -v);
 
