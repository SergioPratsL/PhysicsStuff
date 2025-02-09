function Deriv = DerivadasPotencialYinn( Rx, Ry, v, cad)

% Esta version de la funcion no hace cambio de base, trabaja
% sobre el plano XY con velocidad en direccion X!!
% velocidad normalizada (c=1)

% cad es una cutrada que puede valer t, x, y o z.

    Deriv = [0,0,0,0];
    
    R = sqrt( Rx^2 + Ry^2 );

if cad == 't'
    
    Deriv = [0,0,0,0];
    
elseif cad == 'x'
    
    Deriv(1) = - Rx / R^3;
    Deriv(2) = - Rx / (v * R^3 );
    Deriv(3) = - Ry / (v * R^3);    % 07.12.2015 signo '-' comprobado

    Deriv(4) = 0;
    
elseif cad == 'y'    
    
    Deriv(1) = - Ry / R^3;
    Deriv(2) = - Ry / (v * R^3);
    
% Signo mal!    
%    Deriv(3) = - Rx^3 / (v * R^3 * Ry^2);

% La constante no se deriva respecto la X pero si respecto la Y!!
    dConst = 1 / ( abs(v) *Ry^2);

%    dConst = 9999999;   % Da igual, pase lo que pase, queda inutil...
%    Deriv(3) = Rx^3 / (v * R^3 * Ry^2);    
%    Deriv(3) = dConst + Rx^3 / (v * R^3 * Ry^2);    
    
% 08.12.2015. Uoo!! Debo estar equivocandome ahora pero me da un signo diferente:   
    Deriv(3) = dConst + (Rx/(v*R)) * (1/Ry^2 + 1/R^2);    

    Deriv(4) = 0;    
    
elseif cad == 'z'
    
    Deriv(1) = 0;
    Deriv(2) = 0;
    Deriv(3) = 0;
    
    Constante = - 1 / (abs(v)  * Ry); % Al rotar, la constante sale en la derivada
    Deriv(4) = - Rx / (v * R * Ry^2) + Constante / Ry;  
    
end
