function [Dif_Pot_B, Dif_Pot_A] = VarPotGeodesico_lv2(R_B, v_B, R_C, v_C)

% Package: prueba geodesica nivel 2

% Esta funcion evalua para la particula A lo que cambia el potencila que ve
% de B debido al empujon que le da C en el instante retardado del cual
% estamos recibiendo ahora la luz...

% R_B: vector actual que va de B a A
% v_B: velocidad de B medida desde el SRI de A

% R_C: vector actual que va de C a A
% v_C: velocidad de C medida desde el SRI de A

% Quiero ver la variacion de potencial ocurrida en un intervalo dt medido
% sobre A en el SRI de A.

% Dif_Pot_B: variacion del potencial V_BA por la aceleracion hecha sobre B
% Dif_Pot_A: variacion del potencial V_BA por la aceleracion hecha sobre A


v_0 = [0, 0, 0];
vel_B = norm(v_B);
Sigma_B = 1 / sqrt(1 - vel_B^2);
vel_C = norm(v_C);
Sigma_C = 1 / sqrt(1 - vel_C^2);

% Obtener la posicion retardada de B
Rr_B = PostRetFutEnh( R_B, v_B, 'r');
% Inviertes el sentido del vector porque va de B a A
Rr_B = - Rr_B;

tr = norm(Rr_B);

Rr_B_4 = [tr, Rr_B];

Rr_B_en_B_4 = Boost(Rr_B_4, v_B);

Rr_B_en_B = Rr_B_en_B_4(2:4);       % Distancia recorrida en la posicion antigua


% Boost a la diferencia de posicion entre C y B:
R_CB = - R_C + R_B;     % Estos signos son asi porque quiero un vector de C a B y no al reves

% Vector en el mismo instante con la posicion de C respecto B
R_CB_4 = [0, R_CB];

R_CB_4_en_B = Boost(R_CB_4, v_B);

R_CB_en_B = R_CB_4_en_B(2:4);

% En el SRI de B el los puntos de C y B del "presente" de A no se
% corresponden con el mismo instante en B
t_presente_C = R_CB_4_en_B(1);


% Velocidad de C medida desde B:
v_C_en_B = Vel_Addition_Law( v_C, v_B );

% Averiguar la posicion pasada, ojo con el doble menos!! 
% Necesario porque el vector va de C a B
 Rr_CB_en_B = R_CB_en_B - v_C_en_B * (t_presente_C - tr);
% Rr_CB_en_B = R_CB_en_B - v_C_en_B * (t_presente_C + tr);


% Para saber la posicion desde la que C hizo fuerza a B hay que sacar la
% posicion pasada. Para que vaya bien se ha de poner en sentido de B a C y
% luego dar la vuelta de nuevo
Rold_C_en_B = PostRetFutEnh( - Rr_CB_en_B, v_C_en_B, 'r');
Rold_C_en_B = - Rold_C_en_B;


% Sacar el campo inducido
[E_CB_en_B, B_CB_en_B] = CampoInducido_sin_unidades(Rold_C_en_B, v_C_en_B);

a_B_en_B = FLorentz(v_0, E_CB_en_B, B_CB_en_B);

dv_B_en_B = DifVelRed( a_B_en_B, -v_B, 'E' );

% Desplazamiento espacial para el boost inicial:
dr_B_en_B = - norm(Rr_B_en_B) * a_B_en_B;


% Potencial con el que A ve a B visto desde SRI de B (emisor)
[V_en_B, A_en_B] = PotencialDeYinnAmpliado(Rr_B_en_B, -v_B);

% Primera derivada espacial
[dV_dr, dAx_dr, dAy_dr, dAz_dr] = PrimDerivParcEnh(Rr_B_en_B, -v_B ); 

% Primera derivada de velocidad
[dV_dv, dAx_dv, dAy_dv, dAz_dv] = PrimeraDerivVelEnh(Rr_B_en_B, -v_B);


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


Dif_Pot_n_doppl = Term1 + Term2 + Term3;

 
 
%  ... Aplicar el efecto Doppler
Doppler = Sigma_B * ( 1 - dot( -v_B, Rr_B_en_B / norm(Rr_B_en_B) ) );


% Cuanto mas doppler, mas intenso es el efecto porque significa que toda
% esa aceleracion se recibe en muy poco tiempo.
Dif_Pot_B_en_B = Dif_Pot_n_doppl * Doppler;

Dif_Pot_B = Boost( Dif_Pot_B_en_B, -v_B);


% Tambien hay que sacar la fuerza que C hace sobre B ya que provoca que el
% potencia que B causa en A, V_BA, cambie
Rr_C = PostRetFutEnh( R_C, v_C, 'r');
% Inviertes el sentido del vector porque va de B a A
Rr_C = - Rr_C;
 
tr_C = norm(Rr_C);

Rr_C_4 = [tr_C, Rr_C];

Rr_C_en_B_4 = Boost(Rr_C_4, v_B);

Rr_C_en_B = Rr_C_en_B_4(2:4);   

% Calcular el campo EM que causa C sobre A medido en SRI de B
[E_CA_en_B, B_CA_en_B] = CampoInducido_sin_unidades(Rr_C_en_B, v_C_en_B);

% La dilatacion del tiempo y el aumento de masa inercial se compensan, no
% hay que aplicar ningun factor.
a_A_en_B = FLorentz(-v_B, E_CA_en_B, B_CA_en_B);

dv_A_en_B = DifVelRed( a_A_en_B, -v_B, 'R' );

Dif_Pot_A_en_B(1) = dot(dV_dv, dv_A_en_B);
Dif_Pot_A_en_B(2) = dot(dAx_dv, dv_A_en_B);
Dif_Pot_A_en_B(3) = dot(dAy_dv, dv_A_en_B);
Dif_Pot_A_en_B(4) = dot(dAz_dv, dv_A_en_B);

Dif_Pot_A = Boost(Dif_Pot_A_en_B, -v_B);