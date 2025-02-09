function Deriv2 = SegundDerivPotYinnXX( Rx, Ry, v, dv_e, dv_r, despl )

% Ojo! Para Ay derivo respecto v (aplicando un -(1/v)
% Esto para la constante no es cierto si y=0, pero igualmente eso es algo
% que no probare

R = sqrt( Rx^2 + Ry^2);

Sigma = 1 / sqrt( 1 - v^2);

[ddV, ddAx, ddAy, ddAz] = SecDerivParc(Rx, Ry, v);

[d1_V, d1_Ax, d1_Ay, d1_Az] = PrimDerivParc(Rx, Ry, v);

% [ddV, ddAx, ddAy, ddAz] = DerivEspDerivVel(Rx, Ry, v);

Vect = [Rx, Ry, 0];
[V_yinn, A_yinn] = PotencialDeYinn(Vect, v);



% Efecto(s) del desplazamiento
if despl ~= 0

% V. Efecto del desplazamiento sobre el Minoboost 
    V_pot = V_pot - dv_e * despl * R * ddV(1,1);

% V. Efecto de despl sobre la distancia del Miniboost:
    V_pot = V_pot - dv_e * despl * (Rx/R) * d1_V(1);           
    
% V. Efecto de despl sobre (dV/dv): no hay.

% V. Efecto de despl sobre el boost inverso
    V_pot = V_pot + dv_e * despl * d1_Ax(1);      
        
        
% Ax. Efecto del desplazamiento sobre el Miniboost 
    Ax = Ax - dv_e * despl * R * ddAx(1,1);
      
% Ax. Efecto de despl sobre la distancia del Miniboost:        
    Ax = Ax - dv_e * despl * (Rx/R) * d1_Ax(1);
   
% Ax. Efecto de despl sobre (dV/dv).
    Ax = Ax + dv_e * despl * (1/v) * d1_Ax(1) / Sigma^2;  

% Ax. Efecto de despl sobre el boost inverso
    Ax = Ax + dv_e * despl * d1_V;                  
        
        
% Ay. Efecto desplazamiento sobre el boost inicial
    Ay = Ay - dv_e * despl * R * ddAy(1,1);
     
% Ay. Efecto de despl sobre la distancia del Miniboost:      
    Ay = Ay - dv_e * despl * (Rx/R) * d1_Ay(1);

% Ay. Efecto de despl sobre (dV/dv).    
    Ay = Ay + dv_e * despl * (1/v) * d1_Ay(1) / Sigma^2;  

% Ay. Efecto de despl sobre el boost inverso: no hay (no afecta a Ay)
    
    
% Az: en verdad en este caso vale 0 :)        
     Az = Az + 0;     
    
    
end


% Efecto(s) de la velocidad
if dv_r ~= 0
    
% V. Efecto de dv_r sobre el Miniboost: no hay
        
% V. Efecto de dv_r sobre (dV/dv), no hay.        
        
% V. Efecto de modificar el Sigma que afecta a (dV/dv): no hay

% V. Efecto de dv_r sobre el boost inverso
    V_pot = V_pot - dv_e * dv_r * (1/v) * d1_Ax(1);       

% Ax. Efecto de dv_r sobre el Miniboost
    Ax = Ax + dv_e * dv_r * (R/v)) * d1_Ax(1);      

% Ax. Efecto de dv_r sobre (dAx/dv).
% Derivar velocidad por dv_e pone signo menos y 1/Sigma^2!!
    Ax = Ax - dv_e * dv_r (1/v^2) * A_yinn(1) / Sigma^2; 

% Ax. Efecto de modificar el Sigma que afecta a (dAx/dv): 
    Ax = Ax - dv_e * dv_r * 2 * A_yinn(1);         
        
% Ax. Efecto de dv_r sobre el boost inverso, no hay        
    

% Ay. Efecto de dv_r sobre el Miniboost
    Ay = Ay + dv_e * dv_r * (R/v)) * d1_Ay(1);      

% Ay. Efecto de dv_r sobre (dAx/dv).
    Ay = Ay - dv_e * dv_r (1/v^2) * A_yinn(2) / Sigma^2; 

% Ay. Efecto de modificar el Sigma que afecta a (dAx/dv):     
    Ay = Ay - dv_e * dv_r * 2 * A_yinn(2);      
    
% Ay. Efecto de dv_r sobre el boost inverso: no hay (no afecta a Ay)
     

% Az no cambia     
    Az = Az + 0:

end



Deriv2 = [V_pot, Ax, Ay, Az];