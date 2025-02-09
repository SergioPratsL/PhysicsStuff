function  [V_pot_, Ax_, Ay_, Az_] = SegundDerivPotYinnYY( Rx, Ry, v, dv_e, dv_r, despl )

% 25 terminos :(((

R = sqrt( Rx^2 + Ry^2);

Sigma = 1 / sqrt( 1 - v^2);

[ddV, ddAx, ddAy, ddAz] = SecDerivParc(Rx, Ry, v);

[d1_V, d1_Ax, d1_Ay, d1_Az] = PrimDerivParc(Rx, Ry, v);

% Se necesita para el boost inverso
Vect = [Rx, Ry, 0];
v_vect = [v, 0, 0];
[V_yinn, A_yinn] = PotencialDeYinn(Vect, v_vect);

V_pot = 0;
Ax = 0;
Ay = 0;
Az = 0;

k_rot = 1;


test_1 = ( -Ry * d1_Ax(1) + Rx * d1_Ax(2) )

% Efecto(s) del desplazamiento
if despl ~= 0

% V. Efecto del desplazamiento sobre el Minoboost 
    V_pot = V_pot - dv_e * despl * R * ddV(2,2); 

% V. Efecto de despl sobre la distancia del Miniboost:
    V_pot = V_pot - dv_e * despl * (Ry/R) * d1_V(2);   
    
% V. Efecto de despl sobre (dV/dv)
%    dif_1_V = k_rot * dv_e * despl * (1/v) * ( - Ry*ddV(1,2) + Rx*ddV(1,2) ) / Sigma;   
    dif_1_V = 0;                % 19.02.2016 Derivadas direccionales...   
    dif_2_V = 0;
    V_pot = V_pot + dif_1_V + dif_2_V;

% V. Efecto de despl sobre el boost inverso
    V_pot = V_pot + dv_e * despl * d1_Ay(2);
        
        
% Ax. Efecto del desplazamiento sobre el Miniboost 
    Ax = Ax - dv_e * despl * R * ddAx(2,2); 
      
% Ax. Efecto de despl sobre la distancia del Miniboost:        
    Ax = Ax - dv_e * despl * (Ry/R) * d1_Ax(2);   
    
% Ax. Efecto de despl sobre (dV/dv).
%    dif_1_Ax = k_rot * dv_e * despl * (1/v) * ( - Ry*ddAx(1,2) + Rx*ddAx(1,2) ) / Sigma;   
    dif_1_Ax = 0;               % 19.02.2016 Derivadas direccionales...   
    dif_2_Ax = k_rot * dv_e * despl * (1/v) * d1_Ay(2) / Sigma;
    Ax = Ax + dif_1_Ax + dif_2_Ax;

% Ax. Efecto de despl sobre el boost inverso: no hay     
        
% Ay. Efecto desplazamiento sobre el boost inicial
    Ay = Ay - dv_e * despl * R * ddAy(2,2);
     
% Ay. Efecto de despl sobre la distancia del Miniboost:      
    Ay = Ay - dv_e * despl * (Ry/R) * d1_Ay(2);

% Ay. Efecto de despl sobre (dV/dv).    
    dif_1_Ay = k_rot * dv_e * despl * (1/v) * (- Ry*ddAy(1,2) + Rx*ddAy(1,2) ) / Sigma;   
    dif_2_Ay = - k_rot * dv_e * despl * (1/v) * d1_Ax(2) / Sigma;
    Ay = Ay + dif_1_Ay + dif_2_Ay;

% Ay. Efecto de despl sobre el boost inverso: 
    Ay = Ay + dv_e * despl * d1_V(2);    
    
    
% Az: en verdad en este caso vale 0 :)        
     Az = Az + 0;     
    
end


% Efecto(s) de la velocidad
if dv_r ~= 0
    
% V. Efecto de dv_r sobre el Miniboost
     dif_1_V = k_rot * dv_e * dv_r * (R/v) * ( - Ry * ddV(1,2) + Rx * ddV(2,2) );
    dif_2_V = 0;    % Rotacion entre X e Y no afecta a V.
    V_pot = V_pot + dif_1_V + dif_2_V;

% V. Efecto de dv_r sobre (dV/dv)    (sin k_rot)
%    dif_1_V = - dv_e * dv_r * (1/v^2) * (Rx^2*ddV(2,2) + Ry^2*ddV(1,1) - 2*Rx*Ry*ddV(1,2)) / Sigma;
    div_1_V = 0;                        % 19.02.2016 Derivadas direccionales...    
    dif_2_V = 0;
    V_pot = V_pot + dif_1_V + dif_2_V;

% V. Efecto de modificar el modulo de la velocidad: no hay     
    
% V. Efecto de modificar el Sigma que afecta a (dV/dv): no hay

% V. Efecto de dv_r sobre el boost inverso
    dif_1_V = - k_rot * dv_e * dv_r * (1/v) * ( -Ry*d1_Ay(1) + Rx*d1_Ay(2) );
    dif_2_V = k_rot * dv_e * dv_r * (1/v) * A_yinn(1);
    V_pot = V_pot + dif_1_V + dif_2_V;    

    
% Ax. Efecto de dv_r sobre el Miniboost
    dif_1_Ax = k_rot * dv_e * dv_r * (R/v) * ( -Ry * ddAx(1,2) + Rx * ddAx(2,2) ) ;
    dif_2_Ax = k_rot * dv_e * dv_r * (R/v) * d1_Ay(2);
    Ax = Ax + dif_1_Ax + dif_2_Ax;

% Ax. Efecto de dv_r sobre (dAx/dv).   (sin k_rot el primer termino)
% Fulminando a estos bastardos me acerco al error cero pero deben estar.
%    dif_1_Ax = - dv_e * dv_r * (1/v^2) * ( Rx^2*ddAx(2,2) + Ry^2*ddAx(1,1) - 2*Rx*Ry*ddAx(1,2) ) / Sigma;
    dif_1_Ax = 0;               % 19.02.2016 Derivadas direccionales...       
    dif_2_Ax = + k_rot * 2  * dv_e * dv_r * (1/v^2) * ( -Ry*d1_Ay(1) + Rx*d1_Ay(2)) / Sigma;
    Ax = Ax + dif_1_Ax + dif_2_Ax;

% Ax. Efecto de modificar el modulo de la velocidad  
    Ax = Ax - dv_e * dv_r * (1/v^2) * A_yinn(1) / Sigma^2;    

% Ax. Efecto de modificar el Sigma que afecta a (dAx/dv): no hay 
    
% Ax. Efecto de dv_r sobre el boost inverso: no hay (no afecta a Ax).     
 

% Ay. Efecto de dv_r sobre el Miniboost
    dif_1_Ay = k_rot * dv_e * dv_r * (R/v) * ( -Ry * ddAy(1,2) + Rx * ddAy(2,2) ) ;
    dif_2_Ay = - k_rot * dv_e * dv_r * (R/v) * d1_Ax(2);
    Ay = Ay + dif_1_Ay + dif_2_Ay;    

% Ay. Efecto de dv_r sobre (dAy/dv). (sin k_rot primer termino)
    dif_1_Ay = - dv_e * dv_r * (1/v^2) * ( Rx^2*ddAy(2,2) + Ry^2*ddAy(1,1) - 2*Rx*Ry*ddAy(1,2) ) / Sigma;
    dif_2_Ay = - k_rot * 2 * dv_e * dv_r * (1/v^2) * ( -Ry*d1_Ax(1) + Rx*d1_Ax(2)) / Sigma;
    Ay = Ay + dif_1_Ay + dif_2_Ay;

% Ay. Efecto de modificar el modulo de la velocidad   
    Ay = Ay - dv_e * dv_r * (1/v^2) * A_yinn(2) / Sigma^2;    
    
% Ay. Efecto de modificar el Sigma que afecta a (dAy/dv)   

% Ay. Efecto de dv_r sobre el boost inverso: 
    dif_1_Ay = - k_rot * dv_e * dv_r * (1/v) * ( -Ry*d1_V(1) + Rx*d1_V(2) );
    dif_2_Ay = 0;               % No cambia, V no rota.
    Ay = Ay + dif_1_Ay + dif_2_Ay;    


% Az no cambia     
    Az = Az + 0;

end



%Deriv2 = [V_pot, Ax, Ay, Az];
V_pot_ = V_pot;
Ax_ = Ax;
Ay_ = Ay;
Az_ = Az;