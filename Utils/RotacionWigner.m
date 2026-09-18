function RW = RotacionWigner( v_ab, v_ac );
% Devuelve la matriz de Rotacion resultado de hacer un boost de A a B y
% luego a C.
% Se supone que se parte de A, que conoce v_ab y v_ac, luego aquí
% calculamos v_bc

v_bc = Vel_Addition_Law(v_ac, v_ab);

MB_bc = Tensor_boosts(v_bc);
MB_ab = Tensor_boosts(v_ab);
MB_aci = Tensor_boosts(-v_ac);

RW = MB_bc * MB_ab * MB_aci;