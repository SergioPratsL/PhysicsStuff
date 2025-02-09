function Deriv = SegundDerivPotYinnUnif( Rx, Ry, v, dv_e, dv_r, despl )
% Package: Prueba de la aceleracion dual

% Esta funcion devuelve la variacion de potencial de segudo orden debida a
% la aceleracion del emisor y del receptor, es decir, devuelve el efecto
% de segundo orden provocado por el cambio de veloidad y de posicion del receptor 
% sobre el cambio de potencia que provoca detectar que la velocidad del
% emisor ha cambiado. Es decir saca la diferencia del potencial entre lo
% que habria cambiado el potencial del receptor por a aceleracion del
% emisor en condiciones inerciales del receptor y lo que realmente habra
% cambiado habiando acelerado de forma infinitesimal en cierto momento

% dv_e es la variacion de velocidad del emisor debido a su "acelerin"
% dv_r es la variacion de velocidad del receptor 
% Tanto dv_e como dv_r estan medidas desde el SRI del emisor de forma que
% a dv_r se le ha aplicado la transformacion de velocidades mientras que a

% Ambas son infinitesimales

% v es la velocidad a la que se mueve el receptor, dado que dv_e es
% velocidad propia, se le ha de aplicar transformacion de velocidades para
% sacar la variacion sobre la velocidad v.

% dv_r es la variacion de poscion debida a la variacion de velocidad del
% receptor y al instante en que acelero

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

if dv_e(1) ~= 0
    
% Desplazamiento receptor en X    
    if despl(1) ~= 0 || dv_r(1) ~= 0
        
        [V_pot_aux, Ax_aux, Ay_aux, Az_aux] = SegundDerivPotYinnXX( Rx, Ry, v, dv_e(1), dv_r(1), despl(1) );
        
        V_pot = V_pot + V_pot_aux;
        Ax = Ax + Ax_aux;
        Ay = Ay + Ay_aux;        
        Az = Az + Az_aux;          
        
    end
    
% Desplazamiento receptor en Y   
    if despl(2) ~= 0 || dv_r(2) ~= 0
        
        [V_pot_aux, Ax_aux, Ay_aux, Az_aux] = SegundDerivPotYinnXY( Rx, Ry, v, dv_e(1), dv_r(2), despl(2) );
        
        V_pot = V_pot + V_pot_aux;
        Ax = Ax + Ax_aux;
        Ay = Ay + Ay_aux;        
        Az = Az + Az_aux;          
        
    end
        

end



if dv_e(2) ~= 0

% Desplazamiento receptor en X    
    if despl(1) ~= 0 || dv_r(1) ~= 0
        
        [V_pot_aux, Ax_aux, Ay_aux, Az_aux] = SegundDerivPotYinnYX( Rx, Ry, v, dv_e(2), dv_r(1), despl(1) );
        
        V_pot = V_pot + V_pot_aux;
        Ax = Ax + Ax_aux;
        Ay = Ay + Ay_aux;        
        Az = Az + Az_aux;          
        
    end
    

    if despl(2) ~= 0 || dv_r(2) ~= 0
        
        [V_pot_aux, Ax_aux, Ay_aux, Az_aux] = SegundDerivPotYinnYY( Rx, Ry, v, dv_e(2), dv_r(2), despl(2) );
        
        V_pot = V_pot + V_pot_aux;
        Ax = Ax + Ax_aux;
        Ay = Ay + Ay_aux;        
        Az = Az + Az_aux;          
        
    end
        

end





Deriv = [V_pot, Ax, Ay, Az];