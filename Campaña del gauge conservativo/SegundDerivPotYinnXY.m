function Deriv2 = SegundDerivPotYinnXY( Rx, Ry, v, dv_e, dv_r, despl )




R = sqrt( Rx^2 + Ry^2);

Sigma = 1 / sqrt( 1 - v^2);

[ddV, ddAx, ddAy, ddAz] = SecDerivParc(Rx, Ry, v);

[d1_V, d1_Ax, d1_Ay, d1_Az] = PrimDerivParc(Rx, Ry, v);

% Se necesita para el boost inverso
Vect = [Rx, Ry, 0];
[V_yinn, A_yinn] = PotencialDeYinn(Vect, v);


% Efecto(s) del desplazamiento
if despl ~= 0
    
% V. Efecto del desplazamiento sobre el Minoboost 
    V_pot = V_pot - dv_e * despl * R * ddV(1,2);

% V. Efecto de despl sobre la distancia del Miniboost:
    V_pot = V_pot - dv_e * despl * (Ry/R) * d1_V(1);
    
% V. Efecto de despl sobre (dV/dv)
    V_pot = V_pot + dv_e * despl * (1/v) * d1_V(2) / Sigma^2;

% V. Efecto de despl sobre el boost inverso
    V_pot = V_pot + dv_e * despl * d1_Ax(2);
        
        
% Ax. Efecto del desplazamiento sobre el Miniboost 
    Ax = Ax - dv_e * despl * R * ddAx(1,2);
      
% Ax. Efecto de despl sobre la distancia del Miniboost:        
    Ax = Ax - dv_e * despl * (Ry/R) * d1_Ax(1);
    
% Ax. Efecto de despl sobre (dV/dv).
    Ax = Ax + dv_e * despl * (1/v) * d1_Ax(2) / Sigma^2;

% Ax. Efecto de despl sobre el boost inverso
    Ax = Ax + dv_e * despl * d1_V(2);        


% Ay. Efecto desplazamiento sobre el boost inicial
    Ay = Ay - dv_e * despl * R * ddAy(1,2);
     
% Ay. Efecto de despl sobre la distancia del Miniboost:      
    Ay = Ay - dv_e * despl * (Ry/R) * d1_Ay(1);

% Ay. Efecto de despl sobre (dV/dv).    
    Ay = Ay + dv_e * despl * (1/v) * d1_Ay(2) / Sigma^2;

% Ay. Efecto de despl sobre el boost inverso: no hay (no afecta a Ay)  
    
    
% Az: en verdad en este caso vale 0 :)        
     Az = Az + 0;     
    
end


% Efecto(s) de la velocidad
if dv_r ~= 0
    
% V. Efecto de dv_r sobre el Miniboost
    dif_1_V = - dv_e * dv_r * (1/v) * R * ( -Ry*ddV(1,1) + Rx*ddV(1,2) );
    dif_2_V = 0;
    V_pot = V_pot + dif_1_V + dif_2_V;

% V. Efecto de dv_r sobre (dV/dv)    
    dif_1_V = + dv_e * dv_r * (1/v) * ( -Ry*d1_V(1) + Rx*d1_V(2) ) / Sigma^2;
    dif_2_V = 0;
    V_pot = V_pot + dif_1_V + dif_2_V;    

% V. Efecto de modificar el Sigma que afecta a (dV/dv): no hay


% V. Efecto de dv_r sobre el boost inverso
    dif_1_V = - dv_e * dv_r * (1/v) * ( -Ry*d1_Ax(1) + Rx*d1_Ax(2) );
    dif_2_V = dv_e * dv_r * (1/v) * A_yinn(2);
    V_pot = V_pot + dif_1_V + dif_2_V;        
    
    
% Ax. Efecto de dv_r sobre el Miniboost
    dif_1_Ax = - dv_e * dv_r * (1/v) * R * ( -Ry*ddAx(1,1) + Rx*ddAx(1,2) );
    dif_2_Ax = - dv_e * dv_r * (1/v) * R * d1_Ay(1);
    Ax = Ax + dif_1_Ax + dif_2_Ax;

% Ax. Efecto de dv_r sobre (dAx/dv).
    dif_1_Ax = + dv_e * dv_r * (1/v^2) * ( -Ry*d1_Ax(1) + Rx*d1_Ax(2) ) / Sigma^2;
    dif_2_Ax = + dv_e * dv_r * (1/v^2) * A_yinn(2) / Sigma^2;
    Ax = Ax + dif_1_Ax + dif_2_Ax;    

% Ax. Efecto de modificar el Sigma que afecta a (dAx/dv): no hay 
    

% Ax. Efecto de dv_r sobre el boost inverso: no hay (no afecta a Ax).     
    dif_1_Ax = - dv_e * dv_r * (1/v) * ( -Ry*d1_V(1) + Rx*d1_V(2) );
    dif_2_Ax = 0;
    Ax = Ax + dif_1_Ax + dif_2_Ax;     


% Ay. Efecto de dv_r sobre el Miniboost
    dif_1_Ay = - dv_e * dv_r * (1/v) * R * ( -Ry*ddAy(1,1) + Rx*ddAy(1,2) );
    dif_2_Ay = dv_e * dv_r * (1/v) * R * d1_Ax(1);
    Ay = Ay + dif_1_Ay + dif_2_Ay;

% Ay. Efecto de dv_r sobre (dAy/dv)
    dif_1_Ay = + dv_e * dv_r * (1/v^2) * ( -Ry*d1_Ay(1) + Rx*d1_Ay(2) ) / Sigma^2;
    dif_2_Ay = - dv_e * dv_r * (1/v^2) * A_yinn(1) / Sigma^2;
    Ay = Ay + dif_1_Ay + dif_2_Ay;    


% Ay. Efecto de modificar el Sigma que afecta a (dAy/dv): no hay    


% Ay. Efecto de dv_r sobre el boost inverso: no hay (no afecta a Ay)



% Az no cambia     
    Az = Az + 0;

end



Deriv2 = [V_pot, Ax, Ay, Az];