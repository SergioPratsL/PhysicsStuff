function Deriv = SegundDerivPotYinn( Rx, Ry, v, dv, despl )

% Esta funcion devuelve la variacion de cada potencial en funcion de la
% aceleracion que sufrio el emisor al acelerar y la variacion de posicion
% debida a aceleraciones infinitesimales que el receptor sufrio antes o
% despues que la aceleracion del emisor, pero en todo caso entre tr y tf.
% Es decir, despues de emitir la luz que el emisor recibio en el "0" pero
% antes de recibir la luz que el emisor emitio en el "0"

% dv es la variacion de velocidad del emisor
% despl es el desplazamiento del receptor debido a la aceleracion que tuvo
% el (la del emisor no importa)
% Son cantidades infinitesimales, evaluadas desde el SRI del emisor

% Todo normalizado respecto a la velocidad de la luz

V_pot = 0;      % No cambiara en absoluto
Ax = 0;
Ay = 0;
Az = 0;

R = sqrt( Rx^2 + Ry^2);

if dv(1) ~= 0
    
    if despl(1) ~= 0

        Ax = Ax + dv(1)*despl(1)* Rx/(v^2*R^3);

        Ay = Ay + dv(1)*despl(1)* Ry/(v^2*R^3);
    
        Az = Az + 0;
    
    end
    
    if despl(2) ~= 0    
    
        Ax = Ax + dv(1)*despl(2)* Ry/(v^2*R^3);
        
        Ay = Ay - dv(1)*despl(2)* (Rx / (v^2*R)) * ( (1/Ry^2) + (1/R^2) ) ;
% Termino de la constante:
        Ay = Ay - dv(1)*despl(2)* 1 / (v*abs(v)*Ry^2);
        
        Az = Az + 0;
        
    end        
        
    if despl(3) ~= 0

        Ax = Ax + 0; 

        Ay = Ay + 0; 
    
        Az = Az + dv(1)*despl(3)* Rx / (v^2*R^3);
% Termino de la constante:
        Az = Az + dv(1)*despl(3)* 1 / (v*abs(v)*Ry^2);
        
    end    
    
end    



if dv(2) ~= 0

    if despl(1) ~= 0

        Ax = Ax + dv(2)*despl(1)* Ry/(v^2*R^3);

        Ay = Ay + dv(2)*despl(1)* Rx/(v^2*R^3);
    
        Az = Az + 0;
    
    end
    
    if despl(2) ~= 0

        Ax = Ax - dv(2)*despl(2) * (Rx / (v^2*R)) * ( (1/Ry^2) + (1/R^2) ) ;
% Termino de la constante:
        Ax = Ax - dv(2)*despl(2)* 1 / (v*abs(v)*Ry^2);
        
        Ay = Ay + dv(2)*despl(2)* Ry/(v^2*R^3);
    
        Az = Az + 0;
    
    end    
    
    if despl(3) ~= 0

        Ax = Ax + 0;
        
        Ay = Ay + 0;
    
        Az = Az + 0;    % Curioso caso...
    
    end        
    
end


if dv(3) ~= 0
    
    if despl(1) ~= 0

        Ax = Ax + 0;

        Ay = Ay + 0;
    
        Az = Az + dv(3)*despl(1)* Rx / (v^2*R^3);
    
    end    
    
    if despl(2) ~= 0

        Ax = Ax + 0;

        Ay = Ay + 0;
    
        Az = Az + dv(3)*despl(2)* Ry / (v^2*R^3);
    
    end        
    
    if despl(3) ~= 0

        Ax = Ax + dv(3)*despl(3) * (Rx / (v^2*R)) * ( (1/Ry^2) + (1/R^2) ) ;
% Termino de la constante:
        Ax = Ax + dv(3)*despl(3)* 1 / (v*abs(v)*Ry^2);

        Ay = Ay + 0;
    
        Az = Az + 0;
    
    end            
    
end

Deriv = [V_pot, Ax, Ay, Az];