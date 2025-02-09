function dPot_dv_fut = Deriv_Pot_ret(Rrx, Rry, v, a)

% Esta funcion da la variacion del potencial que "B" ejerce sobre "A"
% debido a una aceleracion de "A".

% El escenario es que estamos en el SRI de "A", " B" se mueve 
% a velocidad v en direccion X. Expresada en relacion a la velocidad de la
% luz.

% Rrx y Rry forma el vector que va de la posicion pasada de "B" a la posicion
% actual de "A" (Al estar en plano XY, Rrz = 0) y obviamente es un vector luz

% a es el vector de aceleracion que sufre "A".

% dPot_dv_fut es la variacion del 4_potencial de "B" debido a la
% aceleracion sufrida por "A"

% v es la velocidad a la que se mueve "B"


% Inicializo:
dif_Pot_dv_x = [0,0,0,0];
dif_Pot_dv_y = [0,0,0,0];
dif_Pot_dv_z = [0,0,0,0];


% Primero de todo obtenemos el Doppler espacial:
Rr = [Rrx, Rry, 0];
v_vect = [v,0,0];

Sigma = 1 / sqrt(1 - v^2);

dist_r = norm(Rr);

% Relacion entre el tiempo que dura un segmento en "A" y lo que tarda "B"
% en observarlo
Doppler = 1 / (1 - v * Rrx / dist_r); 

Rr_t = [dist_r, Rr];

Rr_B_t = Boost(Rr_t, v_vect);

Rr_en_B = [Rr_B_t(2), Rr_B_t(3), Rr_B_t(4)];

% Calculos que aplican a varias componentes
[V_A_desde_B, A_A_desde_B] = PotencialDeYinn(Rr, v_vect);

% 16.12.2015 Ojo! En Z el boost se hace al reves!
% 17.12.2015. Pues todavia seguia mal!
% [V_A_desde_B_paraZ, A_A_desde_B_paraZ] = PotencialDeYinn(Rr, -v_vect);
[V_A_desde_B_paraZ, A_A_desde_B_paraZ] = PotencialDeYinn(Rr_en_B, -v_vect);

% 17.12.2015. Cada componente necesita un potencial diferente...
[V_A_desde_B_paraY, A_A_desde_B_paraY] = PotencialDeYinn(Rr_en_B, -v_vect);


dPot_A_dx = DerivadasPotencialYinn( Rr(1), Rr(2), -v, 'x');

% 16.12.2015. Error... usaba el del script "X"
dPot_A_dx_paraY = DerivadasPotencialYinn( Rr_B_t(2), Rr_B_t(3), -v, 'x');

% 16.12.2015. Error...
%dPot_A_dy = DerivadasPotencialYinn( Rr(1), Rr(2), -v, 'y');
dPot_A_dy = DerivadasPotencialYinn( Rr_B_t(2), Rr_B_t(3), -v, 'y');

% 17.12.2015. Reerror... usar  Rr_B_t(2)
%dPot_A_dz = DerivadasPotencialYinn( Rr(1), Rr(2), -v, 'z');
dPot_A_dz = DerivadasPotencialYinn( Rr_B_t(2), Rr_B_t(3), -v, 'z');


% Trato cada componente por separado.
if a(1) ~= 0    
       
    dV_dvx = 0;
    
    dAx_dvx = - 1 / ((-v)^2 * Rr_B_t(1));    % Rr_b_t(1) es el tiempo, es decir, la distancia    
    
    dAy_dvx = Rr_B_t(2) / ((-v)^2 * Rr_B_t(1) * Rr_B_t(3));    
    dC_dvx =  abs(-v) / ((-v)^3 * Rr_B_t(3));
    dAy_dvx = dAy_dvx + dC_dvx;

    dAz_dvx = 0;
    
    dPot_A_en_B_x = [dV_dvx, dAx_dvx, dAy_dvx, dAz_dvx] / Sigma^2;

    dif_Pot_dv_x = a(1) * Boost( dPot_A_en_B_x,  -v_vect );    
        
end


if a(2) ~= 0

% 16.12.2015 Faltaba el Sigma
    dRot_y = 1 / (v*Sigma);
        
% Rotacion        

    dV_dRot_A_y = - dRot_y * Rr_B_t(3) * dPot_A_dx_paraY(1) + dRot_y * Rr_B_t(2) * dPot_A_dy(1);

    dAx_dRot_A_y = - dRot_y * Rr_B_t(3) * dPot_A_dx_paraY(2) + dRot_y * Rr_B_t(2) * dPot_A_dy(2);

    dAy_dRot_A_y = - dRot_y * Rr_B_t(3) * dPot_A_dx_paraY(3) + dRot_y * Rr_B_t(2) * dPot_A_dy(3);

    dAz_dRot_A_y = - dRot_y * Rr_B_t(3) * dPot_A_dx_paraY(4) + dRot_y * Rr_B_t(2) * dPot_A_dy(4);

    d4V_dRot_A_y = [dV_dRot_A_y, dAx_dRot_A_y, dAy_dRot_A_y, dAz_dRot_A_y];

    
% Desrotamos    

    dV_desrot_A_y = 0;

    dAx_desrot_A_y = dRot_y * A_A_desde_B_paraY(2);

    dAy_desrot_A_y = - dRot_y * A_A_desde_B_paraY(1);

    dAz_desrot_A_y = 0;

    d4V_desrot_A_y = [dV_desrot_A_y, dAx_desrot_A_y, dAy_desrot_A_y, dAz_desrot_A_y];


    dPot_A_desde_B_y = d4V_dRot_A_y + d4V_desrot_A_y;

    
    dif_Pot_dv_y = a(2) * boost(dPot_A_desde_B_y, -v_vect);
    
end


if a(3) ~= 0

    dRot_z = 1 / (v*Sigma);

% Rotacion     
    dV_dRot_A_z = dRot_z * Rr_B_t(2) * dPot_A_dz(1);

    dAx_dRot_A_z = dRot_z * Rr_B_t(2) * dPot_A_dz(2);

    dAy_dRot_A_z = dRot_z * Rr_B_t(2) * dPot_A_dz(3);

    dAz_dRot_A_z = dRot_z * Rr_B_t(2) * dPot_A_dz(4);

    d4V_dRot_A_z = [dV_dRot_A_z, dAx_dRot_A_z, dAy_dRot_A_z, dAz_dRot_A_z];
    
   
% Desrotamos   
    dV_desrot_A_z = 0;

    dAx_desrot_A_z = dRot_z * A_A_desde_B_paraZ(3);

    dAy_desrot_A_z = 0;

    dAz_desrot_A_z = - dRot_z * A_A_desde_B_paraZ(1);

    d4V_desrot_A_z = [dV_desrot_A_z, dAx_desrot_A_z, dAy_desrot_A_z, dAz_desrot_A_z];

    
    dPot_A_desde_B_z = d4V_dRot_A_z + d4V_desrot_A_z;    

    dif_Pot_dv_z = a(3) * boost(dPot_A_desde_B_z, -v_vect);
    
end


% Sumamos las tres componentes y aplicamos el Doppler

% 16.12.2015. Fuera Doppler!!!
% dPot_dv_fut = (dif_Pot_dv_x + dif_Pot_dv_y + dif_Pot_dv_z) / Doppler;
dPot_dv_fut = (dif_Pot_dv_x + dif_Pot_dv_y + dif_Pot_dv_z);