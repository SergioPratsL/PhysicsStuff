function dv_dvr = MatrizConvierteVelocidades(v_vect)
% Package: aceleracion dual


dv_aux = [0,1,0];
dv_aux2 = [0,0,1];

aux = abs( dot( v_vect, dv_aux) );
aux2 = abs( dot( v_vect, dv_aux2) );

if aux < aux2
    dv = dv_aux;
else
    dv = dv_aux2;
end


T = CambioBaseEnh(dv, v_vect)

T_inv = T^(-1);

v = norm(v_vect);
Sigma = 1 / sqrt(1-v^2);

Matriz = [1/Sigma^2, 0, 0; 0, 1/Sigma, 0; 0, 0, 1/Sigma ];

% No se como, pero con v=[0,1,0] necesito esta formula
dv_dvr = T_inv * Matriz * T;



