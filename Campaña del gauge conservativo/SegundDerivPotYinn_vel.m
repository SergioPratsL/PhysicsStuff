function Deriv = SegundDerivPotYinn_vel( Rx, Ry, v, dv_e, dv_r )

% Package: Prueba de la aceleracion dual

% Esta funcion devuelve la variacion de cada potencial en funcion de la
% aceleracion que sufrio el emisor al acelerar y la aceleracion 
% debida a aceleraciones infinitesimales que el receptor sufrio antes o
% despues que la aceleracion del emisor, pero en todo caso entre tr y tf.
% Es decir, despues de emitir la luz que el emisor recibio en el "0" pero
% antes de recibir la luz que el emisor emitio en el "0"

% dv_e es la variacion de velocidad del emisor
% dv_r es la variacion de velocidad del emisor

% Son cantidades infinitesimales, evaluadas desde el SRI del emisor

% Todo normalizado respecto a la velocidad de la luz

V_pot = 0;      
Ax = 0;
Ay = 0;
Az = 0;

R = sqrt( Rx^2 + Ry^2);

Sigma = 1 / sqrt( 1 - v^2);


if dv_e(1) ~= 0
    
    if dv_r(1) ~= 0
        
% 04.01.2016 Segun mis calculos aqui deberia haber un signo menos peor lo descuajeringa todo...
% 15.01.2016. Ojo!! Aqui no va eso!! Esto es un boost inverso sobre el que
% evaluas la variacion de dv_r!!! No invertir el signo aqui!!! Solo si
% ejecutas el (dV/dv_e) !!
        V_pot = V_pot - dv_e(1)*dv_r(1)* 1/(R*v^2);
        
% 04.01.2016 Cambio el signo de esta linea por la "mala nueva"           
        Ax = Ax - dv_e(1)*dv_r(1)* 2/(R*v^3) / Sigma^2;         
        Ax = Ax - dv_e(1)*dv_r(1)* Rx/(R^2*v^2);                
% Termino 09.01.2016
        Ax = Ax - dv_e(1)*dv_r(1)* 2/(v*R);
        
% Golpe criminal para hacer que una prueba salga bien... meter el cuadrado
% en la redonda a base de hostias :P.
%        Ax = Ax - abs(v) / v * 9.925400698
        
        
% 04.01.2016 Cambio el signo de esta linea por la "mala nueva"    
        Ay = Ay - dv_e(1)*dv_r(1)* ( 2*Rx/(R*Ry*v^3) + 2/(Ry*v^2*abs(v)) ) / Sigma^2;
        Ay = Ay - dv_e(1)*dv_r(1)* Ry/(R^2*v^2);
% Termino 09.01.2016
        Ay = Ay + dv_e(1)*dv_r(1)* (2/Ry)*( Rx/(v*R) + 1/abs(v)  );
        
        Az = Az + 0;
        
    end
    
end


% En construccion!!

Deriv = [V_pot, Ax, Ay, Az];