function T_boost = BoostTensorEE(T, u)
% Lambda es las matriz de transformación de los boosts

    Lambda = Tensor_boosts(u);
    %Lambda_i = Tensor_boosts(-u);

    % ¡Fuerza bruta al ataque!
    %EE_boost = Lambda * Lambda * EE;
    % Lambda = Lambda' --> perdida de tiempo
    %EE_boost = Lambda' * Lambda * EE;  
    %EE_boost = Lambda * EE;    % Ni hablar!
    %EE_boost = Lambda * EE * Lambda_i;
    %EE_boost = Lambda_i * EE * Lambda;
    
    % A veces la fuerza bruta es el camino mas rapido hacia la "recta
    % opinion"
    T_boost = Lambda * T * Lambda;
end