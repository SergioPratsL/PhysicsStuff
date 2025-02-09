function dPot_dv_fut = Deriv_Pot_fut(Rfx, Rfy, v, a)

% Esta funcion da la variacion del potencial que "A" ejerce sobre "B"
% debido a una aceleracion de "A".

% El escenario es que estamos en el SRI de "A", " B" se mueve 
% a velocidad v en direccion X. Expresada en relacion a la velocidad de la
% luz.

% Rfx y Rfy forma el vector que va de la posicion actual de "A" a la posicion
% futura de B (Al estar en plano XY, Rfz = 0) y obviamente es un vector luz

% a es el vector de aceleracion que sufre "A".

% dPot_dv_fut es la variacion del 4_potencial de "B" debido a la
% aceleracion sufrida por "A"

% v es la velocidad a la que se mueve "B"

% Inicializo:
dif_Pot_dv_x = [0,0,0,0];
dif_Pot_dv_y = [0,0,0,0];
dif_Pot_dv_z = [0,0,0,0];


% Primero de todo obtenemos el Doppler espacial:
Rf = [Rfx, Rfy, 0];
v_vect = [v,0,0];

Sigma = 1 / sqrt(1 - v^2);

dist_f = norm(Rf);

% Relacion entre el tiempo que dura un segmento en "A" y lo que tarda "B"
% en observarlo
Doppler = 1 * (1 - v * Rfx / dist_f);   % OK


% Calculos que aplican a varias componentes
[V_B_desde_A, A_B_desde_A] = PotencialDeYinn(Rf, v_vect);

dPot_B_dx = DerivadasPotencialYinn( Rf(1), Rf(2), v, 'x');

dPot_B_dy = DerivadasPotencialYinn( Rf(1), Rf(2), v, 'y');

dPot_B_dz = DerivadasPotencialYinn( Rf(1), Rf(2), v, 'z');


% Trato cada componente por separado.
if a(1) ~= 0    
    
    dV_dvx = 0;
    dAx_dvx = - 1 / ( v^2 * dist_f );
    dAy_dvx = Rf(1) / ( v^2 * dist_f * Rf(2) );

    dC_dvx = abs(v) / ( v^3 * Rf(2));
    dAy_dvx = dAy_dvx + dC_dvx;

    dAz_dvx = 0;
    
% Si acelero el origen la aceleracion es en sentido contrario (signo
% menos)! Tambien debe aplicars el factor 1/Sigma^2
    dPot_B_dvx = - [dV_dvx, dAx_dvx, dAy_dvx, dAz_dvx] / Sigma^2;

    dPot_Boost_final_x =  [A_B_desde_A(1), V_B_desde_A, 0, 0];
    
    dif_Pot_dv_x = a(1) * ( - dist_f * dPot_B_dx + dPot_B_dvx + dPot_Boost_final_x);
    
end


if a(2) ~= 0

% Boost inicial 
    dPot_B_mini_boost_y = - dist_f * dPot_B_dy;
    
    
% Rotacion    
    dRot_y = 1 / v;
    
    dV_dRot_B_y = - dRot_y * Rf(2) * dPot_B_dx(1) + dRot_y * Rf(1) * dPot_B_dy(1);

    dAx_dRot_B_y = - dRot_y * Rf(2) * dPot_B_dx(2) + dRot_y * Rf(1) * dPot_B_dy(2);

    dAy_dRot_B_y = - dRot_y * Rf(2) * dPot_B_dx(3) + dRot_y * Rf(1) * dPot_B_dy(3);   % Parece OK

    dAz_dRot_B_y = - dRot_y * Rf(2) * dPot_B_dx(4) + dRot_y * Rf(1) * dPot_B_dy(4);

    d4V_dRot_B_y = [dV_dRot_B_y, dAx_dRot_B_y, dAy_dRot_B_y, dAz_dRot_B_y];    

    
% Desrotamos     
    dV_desrot_B_y = 0;

    dAx_desrot_B_y = dRot_y * A_B_desde_A(2);

    dAy_desrot_B_y = - dRot_y * A_B_desde_A(1);

    dAz_desrot_B_y = 0;

    d4V_desrot_B_y = [dV_desrot_B_y, dAx_desrot_B_y, dAy_desrot_B_y, dAz_desrot_B_y];

    
% Deshacemos el mini-boost    
    dV_boost_inv_B_y = A_B_desde_A(2);

    dAx_boost_inv_B_y = 0;

    dAy_boost_inv_B_y = V_B_desde_A;

    dAz_boost_inv_B_y = 0;

    dPot_B_boost_inv_y = [dV_boost_inv_B_y, dAx_boost_inv_B_y, dAy_boost_inv_B_y, dAz_boost_inv_B_y];
    
    
    dif_Pot_dv_y = a(2) * (dPot_B_mini_boost_y + d4V_dRot_B_y + d4V_desrot_B_y + dPot_B_boost_inv_y);
    
end


if a(3) ~= 0

% Boost inicial   
    dPot_B_mini_boost_z = - dist_f * dPot_B_dz;


% Rotacion     
    dRot_z = 1 / v;    
 
    dV_dRot_B_z = dRot_z * Rf(1) * dPot_B_dz(1);

    dAx_dRot_B_z = dRot_z * Rf(1) * dPot_B_dz(2);

    dAy_dRot_B_z = dRot_z * Rf(1) * dPot_B_dz(3);

    dAz_dRot_B_z = dRot_z * Rf(1) * dPot_B_dz(4);

    d4V_dRot_B_z = [dV_dRot_B_z, dAx_dRot_B_z, dAy_dRot_B_z, dAz_dRot_B_z];    
    
    
% Desrotamos   
    dV_desrot_B_z = 0;

    dAx_desrot_B_z = dRot_z * A_B_desde_A(3);

    dAy_desrot_B_z = 0;

    dAz_desrot_B_z = - dRot_z * A_B_desde_A(1);

    d4V_desrot_B_z = [dV_desrot_B_z, dAx_desrot_B_z, dAy_desrot_B_z, dAz_desrot_B_z];
    
    
% Deshacemos el mini-boost
    dV_boost_inv_B_z = A_B_desde_A(3);  % Entiendo que sera 0...

    dAx_boost_inv_B_z = 0;

    dAy_boost_inv_B_z = 0;

    dAz_boost_inv_B_z = V_B_desde_A;

    dPot_B_boost_inv_z = [dV_boost_inv_B_z, dAx_boost_inv_B_z, dAy_boost_inv_B_z, dAz_boost_inv_B_z];
    
    
    dif_Pot_dv_z = a(3) * (dPot_B_mini_boost_z + d4V_dRot_B_z + d4V_desrot_B_z + dPot_B_boost_inv_z);
    
end


% Sumamos las tres componentes y aplicamos el Doppler

% 16.12.2015. Fuera Doppler!!!
%dPot_dv_fut = (dif_Pot_dv_x + dif_Pot_dv_y + dif_Pot_dv_z) * Doppler;
dPot_dv_fut = (dif_Pot_dv_x + dif_Pot_dv_y + dif_Pot_dv_z);
