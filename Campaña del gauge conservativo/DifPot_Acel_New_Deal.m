function [Pot_A, Pot_B, F_rad, Dif] = DifPot_Acel_New_Deal( R_vect, v_vect, a)
% R_vect representa la posicion actual de B respecto a A
% v_vect representa la velocidad de B respecto el SRI de A

% El que acelera es A

% La aceleracion de ahora ve como se modifica el potencial procedente de la
% posicion pasada de B, por otro lado B notara el efecto de la aceleracion
% en el instante futuro sujeto al retardo de propagacion de la luz que sale
% ahora

R = norm(R_vect);
v = norm(v_vect);

Sigma = 1 / sqrt(1 - v^2);



%Rr_vect = PostRetFut( R_vect, v_vect(1), 'r');
Rr_vect = PostRetFutEnh(R_vect, v_vect, 'r');
% Inviertes el sentido del vector
Rr_vect = - Rr_vect;

% Rf_vect = PostRetFut( R_vect, v_vect(1), 'f');
Rf_vect = PostRetFutEnh(R_vect, v_vect, 'f');

Rr = norm(Rr_vect);
Rf = norm(Rf_vect);

Rr_vect_t = [Rr, Rr_vect];
Rr_vect_B_t = Boost(Rr_vect_t, v_vect);

Rr_vect_en_B = Rr_vect_B_t(2:4);


% Variacion del potencial que ve A, debe ser evaluada desde B y luego hacer
% boost, por ello solo depende de la variacion dv

% Obtenemos la aceleracion vista desde el SRI de B:
a_en_B = DifVelRed( a, -v_vect, 'R' );


% Calculo formalmente las derivadas
[dV_dv_A, dAx_dv_A, dAy_dv_A, dAz_dv_A] = PrimeraDerivVelEnh(Rr_vect_en_B, -v_vect);

dPot_A_en_B(1) = dot(dV_dv_A, a_en_B);
dPot_A_en_B(2) = dot(dAx_dv_A, a_en_B);
dPot_A_en_B(3) = dot(dAy_dv_A, a_en_B);
dPot_A_en_B(4) = dot(dAz_dv_A, a_en_B);

Pot_A = Boost( dPot_A_en_B,  -v_vect );



% Variacion del potencial que ve B (en el futuro) por la aceleracion de A
% Se evalua desde el SRI de A pero requiere tres terminos

[dV_dr_B, dAx_dr_B, dAy_dr_B, dAz_dr_B] = PrimDerivParcEnh(Rf_vect, v_vect); 

[dV_dv_B, dAx_dv_B, dAy_dv_B, dAz_dv_B] = PrimeraDerivVelEnh(Rf_vect, v_vect); 

[V_B_aux, A_B_aux] = PotencialDeYinn(Rf_vect, v_vect);       % Solo para comparar
[V_B, A_B] = PotencialDeYinnAmpliado(Rf_vect, v_vect);


Term_1(1) = - Rf * dot(dV_dr_B, a);
Term_1(2) = - Rf * dot(dAx_dr_B, a);
Term_1(3) = - Rf * dot(dAy_dr_B, a);
Term_1(4) = - Rf * dot(dAz_dr_B, a);

dif_v = DifVelRed( a, -v_vect, 'E' );

Term_2(1) = dot(dV_dv_B, dif_v);
Term_2(2) = dot(dAx_dv_B, dif_v);
Term_2(3) = dot(dAy_dv_B, dif_v);
Term_2(4) = dot(dAz_dv_B, dif_v);

Term_3(1) = dot(a, A_B);
Term_3(2) = V_B * a(1);
Term_3(3) = V_B * a(2);
Term_3(4) = V_B * a(3);


Pot_B = Term_1 + Term_2 + Term_3;


v_0 = [0, 0, 0];    % Se evalua desde el emisor.

[E, B] = CampoRadiadoCargaAcelerada(Rf_vect, v_0, a);

Dopler = 1 - dot(v_vect, Rf_vect) / Rf;

F = E + cross(v_vect, B);


% Warning!!
F_rad = [dot(F, v_vect), F] / Dopler;


Dif = Pot_A + Pot_B + F_rad;