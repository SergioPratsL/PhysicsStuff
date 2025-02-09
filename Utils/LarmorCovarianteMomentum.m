function Pot_F = LarmorCovarianteMomentum( v, acel, Sigma )
% Esto obtiene la potencia radiada, basado en fórmula de Richard
% Fiztpatrick
% https://farside.ph.utexas.edu/teaching/em/lectures/node130.html
% v y acel son la velocidad y aceleración medidas desde el sistema de
% reposo, siempre con la velocidad de la luz igual a 1.
% El 4piEpsilon también igual a 1, la carga también igual a 1. Make it easy
% Si sigma ya está calculada por la función que llama, que la pase, si no,
% que pase un 0 y se calcula aquí.

    if( Sigma == 0)
        Sigma = fGamma(norm(v));
    end

    Sigma6 = Sigma^6;
    
    Pot = (2/3) * Sigma6 * (norm(acel)^2 - norm(cross(v, acel))^2);
    
    % Por ser un boost de un vector original de tipo energia [1,0,0,0].
    % [1,0,0,0] --> Sigma * [1, v]
    F = v * Pot;    % Correcta
    %F = v * Pot * Sigma;   % BAD
    
    Pot_F = [Pot, F];
    
end

