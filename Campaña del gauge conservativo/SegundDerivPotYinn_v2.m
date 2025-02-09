function Deriv = SegundDerivPotYinn_v2( Rx, Ry, v, dv_e, despl )

% Package: Prueba de la aceleracion dual

% Esta funcion devuelve la variacion de cada potencial en funcion de la
% aceleracion que sufrio el emisor al acelerar y la variacion de posicion
% debida a aceleraciones infinitesimales que el receptor sufrio antes o
% despues que la aceleracion del emisor, pero en todo caso entre tr y tf.
% Es decir, despues de emitir la luz que el emisor recibio en el "0" pero
% antes de recibir la luz que el emisor emitio en el "0"

% dv_e es la variacion de velocidad del emisor
% despl es el desplazamiento del receptor debido a la aceleracion que tuvo
% (la del emisor no importa)

% Son cantidades infinitesimales, evaluadas desde el SRI del emisor

% Todo normalizado respecto a la velocidad de la luz

V_pot = 0;     
Ax = 0;
Ay = 0;
Az = 0;

% Sacare todas las derivadas espaciales... no me servira para las segunda
% derivadas...
DPot_x = DerivadasPotencialYinn( Rx, Ry, v, 'x');

DPot_y = DerivadasPotencialYinn( Rx, Ry, v, 'y');

DPot_z = DerivadasPotencialYinn( Rx, Ry, v, 'z');


R = sqrt( Rx^2 + Ry^2);

Sigma = 1 / sqrt( 1 - v^2);
% Sigma = 1; Con esto el error aumenta :)

if dv_e(1) ~= 0
    
    if despl(1) ~= 0

% 15.01.2015. ¿Como puse este termino???
        V_pot = V_pot + dv_e(1)*despl(1)* R / (Ry^2/R^5 - 2*Rx^2/R^5);
% Nuevo 08.01.2016 (pues esta mal pero no lo tocare)
        V_pot = V_pot + dv_e(1)*despl(1) * Rx^2/R^4;
        V_pot = V_pot + dv_e(1)*despl(1)*DPot_x(2);  % dAx/dx
        
% 04.01.2016 Cambio el signo de esta linea por la "mala nueva"        
        Ax = Ax - dv_e(1)*despl(1)* Rx/(v^2*R^3) / Sigma^2;
% Tb cambia de signo por "error humano"        
        Ax = Ax + dv_e(1)*despl(1)* (R/v) * (Ry^2/R^5 - 2*Rx^2/R^5);
% Nuevo 08.01.2016.        
        Ax = Ax + dv_e(1)*despl(1) * Rx^2/(v*R^4);
% Este no ha de cambiar de signo, DPot_x es de fiar :)
        Ax = Ax + dv_e(1)*despl(1)* DPot_x(1);      % dV/dx

% 04.01.2016 Cambio el signo de esta linea por la "mala nueva" 
% 04.01.2016. Otro fallo garrafal: Ry/(v^2*R^3)#^2#  (sobra el ^2)
        Ay = Ay - dv_e(1)*despl(1)* Ry/(v^2*R^3) / Sigma^2;      
        Ay = Ay - dv_e(1)*despl(1)* (R/v) * (3*Rx*Ry/R^5);      % OK
    
        Az = Az + 0;
    
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    if despl(2) ~= 0    
    
        V_pot = V_pot - dv_e(1)*despl(2) * R * (3*Rx*Ry/R^5);
        V_pot = V_pot + dv_e(1)*despl(2)*DPot_y(2);  % dAx/dy        
        
        Ax = Ax + dv_e(1)*despl(2)* Ry/(v^2*R^3);
        Ax = Ax - dv_e(1)*despl(2) * (R/v) * (3*Rx*Ry/R^5);
        Ax = Ax + dv_e(1)*despl(2) * DPot_y(1);      % dV/dy
        
        Ay = Ay - dv_e(1)*despl(2)* (Rx / (v^2*R)) * ( (1/Ry^2) + (1/R^2) ) ;
        Ay = Ay - dv_e(1)*despl(2)* 1/(v*R) * (1/Ry^2 + Ry^2/R^4 - Rx^2/R^4);
% Termino de la constante:
        Ay = Ay - dv_e(1)*despl(2)* 1 / (v*abs(v)*Ry^2);
        
        Az = Az + 0;
        
    end        
        
    if despl(3) ~= 0
        
        V_pot = V_pot + 0;

        Ax = Ax + 0; 

        Ay = Ay + 0; 
    
        Az = Az + dv_e(1)*despl(3)* Rx / (v^2*R^3);
% Termino de la constante:
        Az = Az + dv_e(1)*despl(3)* 1 / (v*abs(v)*Ry^2);
        
    end    
    
end    



if dv_e(2) ~= 0

    if despl(1) ~= 0

        V_pot = V_pot - dv_e(2)*despl(1)* R* (3*Rx*Ry/R^5);        
        V_pot = V_pot + dv_e(2)*despl(1)* DPot_x(3);  % dAy/dx
        
        Ax = Ax + dv_e(2)*despl(1) * Ry/(v^2*R^3);
        Ax = Ax - dv_e(2)*despl(1) * (R/v) * (3*Rx*Ry/R^5);

        Ay = Ay + dv_e(2)*despl(1) * Rx/(v^2*R^3);
        Ay = Ay - dv_e(2)*despl(1) * (R/v) * (Rx^2/R^5 - 2*Ry^2/R^5);
        Ay = Ay + dv_e(2)*despl(1) * DPot_x(1);         % dV/dx
        
        Az = Az + 0;
    
    end
    
    if despl(2) ~= 0

        V_pot = V_pot - dv_e(2)*despl(2) * R * (Rx^2/R^5 - 2*Ry^2/R^5);
        
        V_pot = V_pot + dv_e(2)*despl(2) * DPot_y(3);  % dAy/dy
              
        Ax = Ax - dv_e(2)*despl(2) * (Rx / (v^2*R)) * ( (1/Ry^2) + (1/R^2) ) ;
% Termino de la constante:
        Ax = Ax - dv_e(2)*despl(2)* 1 / (v*abs(v)*Ry^2);
        Ax = Ax - dv_e(2)*despl(2) - (R/v) * (Rx^2/R^5 - 2*Ry^2/R^5);
        
        Ay = Ay + dv_e(2)*despl(2)* Ry/(v^2*R^3);
        Ay = Ay - dv_e(2)*despl(2)* R * (-2)*( 1/(abs(v)*R) + (Rx/v*R) * (1/Ry^3 + Ry/R^4) );
        Ay = Ay + dv_e(2)*despl(2)* DPot_y(1);         % dV/dy
    
        Az = Az + 0;
    
    end    
    
    if despl(3) ~= 0
        
        V_pot = V_pot + 0;

        Ax = Ax + 0;
        
        Ay = Ay + 0;
    
        Az = Az + 0;    % Curioso caso...
    
    end        
    
end


if dv_e(3) ~= 0
    
% En construccion!
    auxy = 1 / 0;
    
    if despl(1) ~= 0

        V_pot = V_pot + 0;
        
        Ax = Ax + 0;

        Ay = Ay + 0;
    
        Az = Az + dv_e(3)*despl(1)* Rx / (v^2*R^3);
    
    end    
    
    if despl(2) ~= 0

        V_pot = V_pot + 0;        
        
        Ax = Ax + 0;

        Ay = Ay + 0;
    
        Az = Az + dv_e(3)*despl(2)* Ry / (v^2*R^3);
    
    end        
    
    if despl(3) ~= 0

        V_pot = V_pot + dv_e(3)*despl(3) * (1/R^2 +  DPot_z(4) );  % dAz/dz
        
        Ax = Ax + dv_e(3)*despl(3) * (Rx / (v^2*R)) * ( (1/Ry^2) + (1/R^2) ) ;
% Termino de la constante:
        Ax = Ax + dv_e(3)*despl(3)* 1 / (v*abs(v)*Ry^2);
        Ax = Ax - dv_e(3)*despl(3)* (-1)/(R^2*v);
        
        Ay = Ay - dv_e(3)*despl(3) * Rx/(v*Ry) * 1/R^2;
        
    
        Az = Az + 0;
    
    end            
    
end

Deriv = [V_pot, Ax, Ay, Az];