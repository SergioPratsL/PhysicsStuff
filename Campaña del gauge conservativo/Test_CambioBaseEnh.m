
% Pruena 1. OK.
R_vect = [1; 2; 0];
v_vect = [1; 0; 0];


% Prueba 2. OK
R_vect = [1; 2; -1];
v_vect = [1; 1; 0];


% Prueba 3. OK
R_vect = [4; 2; -1];
v_vect = [0; 0.54; 0.19];


% Suficiente :).


T = CambioBaseEnh(R_vect, v_vect)

T_inv = T^(-1)

R_rot = T * R_vect

v_rot = T * v_vect

R_original =  T_inv * R_rot

v_original =  T_inv * v_rot


