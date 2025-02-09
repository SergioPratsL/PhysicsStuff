
% Prueba piloto, esta se supera finalmente. OK
% R_vect = [0, 1, 0];
% v_vect = [0.4, 0, 0];
% a = [1, 0, 0];

% Prueba 2
% R_vect = [1, 2, 0];
% v_vect = [-0.6, 0, 0];
% a = [1, 0, 0];


% Prueba 3
% R_vect = [0, 1, 0];
% v_vect = [0.4, 0, 0];
% a = [0, 1, 0];


% Prueba 4. OK (la primera que no ha requerido modificaciones
% R_vect = [2, -1, 0];
% v_vect = [0.4, 0, 0];
% a = [0, -1, 0];


% Prueba 5. OK
% R_vect = [2, -1, 0];
% v_vect = [0.3, 0, 0];
% a = [0.5, -1, 0];


% Prueba 6. OK, :D
% R_vect = [2, 1.5, 0];
% v_vect = [-0.3, 0, 0];
% a = [0, -0.4, 1];


% Prueba 7. OK
% R_vect = [-0.5, 0.4, 0.9];
% v_vect = [0, 0.4, 0];
% a = [0.2, -0.4, 0.6];


% Prueba 8. 
R_vect = [-0.5, 0.4, 0.9];
v_vect = [0, 0, 0.4];
a = [0.2, -0.4, 0.6];


% Prueba 9. 
%  R_vect = [-0.5, 0.4, 0.9];
%  v_vect = [0.4, 0, 0];
%  a = [0.2, -0.4, 0.6];




[Pot_A, Pot_B, F_rad, Dif] = DifPot_Acel_New_Deal( R_vect, v_vect, a)