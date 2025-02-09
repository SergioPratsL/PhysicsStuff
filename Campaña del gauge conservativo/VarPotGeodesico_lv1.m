function [Dif_Pot_A, Dif_Pot_B, Rad] = VarPotGeodesico_lv1(R_vect, v_vect)

% Package: prueba geodesica nivel 1

% Esta funcion evalua para la particula A lo que 
% cambia el potencial en cierto instante debido a le 
% interaccion con otra particula B teniendo en cuenta tres factores:

% Como cambia el potencial por la aceleracion de A
% Como cambia el potencial por la aceleracion de B que se debe a la fuerza
% retardada de A
% Fuerza y potencia (potencia es 0) causada por el campo radiado de B

% Tratamos con particulas infinitesimales con lo que todos los terminos son
% neglibilbes pero quiero evaluar si las variaciones van parejas
% La carga y la masa son dV^3 y para simplificar asumo que dq / dm = 1 es
% decir establezco la unidad entre carga y masa.

% El potencia es dV^6 asi como la fuerza del campo inducido.

% Comparamos un dt de recepcion en A con el tiempo asociado en B (afectado
% por la inversa del efecto Doppler

% Vaya lio!! Se tiene que evaluar todo desde el SRI de B y solo al final
% pasar al SRI de A.

% R_vect va de B a A, al contrario de lo que hice en los test de
% aceleracion originales (lo he liado todo esta vez, B ocupa en verdad el
% lugar de A y viceversa).


v_0 = [0, 0, 0];
R = norm(R_vect);
v = norm(v_vect);
Sigma = 1 / sqrt(1 - v^2);

R_vect_4 = [R, R_vect];

R_vect_en_B_4 = Boost(R_vect_4, v_vect);

% R_vect_en_B sigue yendo de B a A pero es visto desde el SRI de B
R_vect_en_B = R_vect_en_B_4(2:4);
R_en_B = norm(R_vect_en_B);


% Saco la variacion de la velocidad con este pesado pero seguro metodo:
% [dS_dvi, dS_dvr] = PrimDerivLeyAddVeloc(-v_vect, v_0);
% Chufa de dS_dvr, asume que dvr es medido desde SRI de A!!



% Potencial con el que A ve a B visto desde SRI de B (emisor)
[V_en_B, A_en_B] = PotencialDeYinnAmpliado(R_vect_en_B, -v_vect);


% Primera derivada espacial
[dV_dr, dAx_dr, dAy_dr, dAz_dr] = PrimDerivParcEnh(R_vect_en_B, -v_vect); 

% Primera derivada de velocidad
[dV_dv, dAx_dv, dAy_dv, dAz_dv] = PrimeraDerivVelEnh(R_vect_en_B, -v_vect);



% Doppler inverso, necesario para relacionar los intervalos del emisor y el
% receptor
% B es realmente quien emite y cogemos un intervalo dt', el intervalo del
% receptor crece a razon del doppler inverso
 Dopler_inv_en_B = 1 / (1 - dot( R_vect_en_B / R_en_B, - v_vect)  );


% Vemos el campo que genera B sobre A (en SRI de A)
 [E_BA, B_BA] = CampoInducido_sin_unidades(R_vect, v_vect);

 % Campo que genera B sobre A visto desde el SRI de B
 [E_BA_en_B, B_BA_en_B] = Boost_EM(E_BA, B_BA, v_vect);
 
% Deberia valor lo mismo que E_BA_en_B
 aux = R_vect_en_B / R_en_B^3;      % OK!
 
 % El sigma es por el aumento de la masa inercial
 a_A_en_B = FLorentz(-v_vect, E_BA_en_B, B_BA_en_B) / Sigma;
 
 dv_A_aux = a_A_en_B * Dopler_inv_en_B;

 dv_A_en_B = DifVelRed( dv_A_aux, -v_vect, 'R' );
 


% Vemos el campo que genera A sobre B (en el SRI de A)
[E_AB, B_AB] = CampoInducido_sin_unidades(-R_vect, v_0);

% Campo que genera A sobre B visto desde el SRI de B
 [E_AB_en_B, B_AB_en_B] = Boost_EM(E_AB, B_AB, v_vect);
 

 
% Aceleracion de B:
a_B_en_B = FLorentz(v_0, E_AB_en_B, B_AB_en_B);

% dt', el tiempo de B, pasa a ser el unitario
dv_B_aux = a_B_en_B;
dv_B_en_B = DifVelRed( dv_B_aux, -v_vect, 'E' );

% Necesitamos el vector que va de B a A en el SRI de B
% Ya lo tenemos, es R_vect_en_B...
dr_B_en_B = - R_en_B * a_B_en_B;



% % Sacare el campo visualizado en el SRI de B para comparar.
% [E_AB_aux, B_AB_aux] = Boost_EM(E_AB, B_AB, v_vect);
% 
% % [kk1, kk2] = Boost_EM(E_AB_aux, B_AB_aux, -v_vect); % OK
% 
% a_B_aux = FLorentz(v_0, E_AB_aux, B_AB_aux);
% 
% % Divido por Sigma porque el tiempo propio es Sigma veces mas corto
% dV_B_aux = (dS_dvi * a_B_aux') * Dopler_inv / Sigma; 
% 
% % Otro camino
% vect_aux = [R, -R_vect];
% vect_aux_transf = Boost(vect_aux, v_vect);
% R_aux = [vect_aux_transf(2), vect_aux_transf(3), vect_aux_transf(4)];
% [E_AB_aux2, B_AB_aux2] = CampoInducido_sin_unidades(R_aux, -v_vect);
% 
% a_B_aux2 = FLorentz(v_0, E_AB, B_AB) ; 
% 
% dV_B_aux2 = (dS_dvi * a_B_aux2') * Dopler_inv / Sigma;  % Sigma por reducir el tiempo


% Campo radiado visto desde B, como no
[E_rad_en_B, B_rad_en_B] = CampoRadiadoCargaAcelerada( R_vect_en_B, v_0, dv_B_en_B);
F_rad_en_B = FLorentz(v_0, E_rad_en_B, B_rad_en_B);

% El efecto de la radiacion se alarga Doppler_inv_en_B!!
Rad_en_B = [dot(-v_vect, F_rad_en_B), F_rad_en_B] * Dopler_inv_en_B;


Dif_Pot_A_en_B(1) = dot(dV_dv, dv_A_en_B);
Dif_Pot_A_en_B(2) = dot(dAx_dv, dv_A_en_B);
Dif_Pot_A_en_B(3) = dot(dAy_dv, dv_A_en_B);
Dif_Pot_A_en_B(4) = dot(dAz_dv, dv_A_en_B);


% La diferencia por la aceleracion del emisor tiene tres terminos:
% variacion de posicion (por boost) -> Term1
% variacion de velocidad            -> Term2
% boost inverso                     -> Term3
Term1(1) = dot(dV_dr, dr_B_en_B);
Term1(2) = dot(dAx_dr, dr_B_en_B);
Term1(3) = dot(dAy_dr, dr_B_en_B);
Term1(4) = dot(dAz_dr, dr_B_en_B);

Term2(1) = dot(dV_dv, dv_B_en_B);
Term2(2) = dot(dAx_dv, dv_B_en_B);
Term2(3) = dot(dAy_dv, dv_B_en_B);
Term2(4) = dot(dAz_dv, dv_B_en_B);

Term3(1) = dot(A_en_B, a_B_en_B);
Term3(2) = V_en_B * a_B_en_B(1);
Term3(3) = V_en_B * a_B_en_B(2);
Term3(4) = V_en_B * a_B_en_B(3);


Dif_Pot_B_en_B = Term1 + Term2 + Term3;



Dif_Pot_A = Boost( Dif_Pot_A_en_B, -v_vect);

Dif_Pot_B = Boost( Dif_Pot_B_en_B, -v_vect);

Rad = Boost( Rad_en_B, -v_vect);
