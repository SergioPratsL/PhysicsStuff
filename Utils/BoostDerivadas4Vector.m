function D_A_boost = BoostDerivadas4Vector(D_A, u)
% Lambda es las matriz de transformación de los boosts
    Lambda = Tensor_boosts(u);
    Lambda_i = Tensor_boosts(-u);

    D_A_boost = Lambda * D_A * Lambda_i;
    % Parece "conveniente", pero no es el correcto.
    %D_A_boost = Lambda_i * D_A * Lambda;
end


% Pruebas para saber si BoostDerivadas4Vector es correcto, parece que si
%  D_A = [0,0,0,0; 0,0,0,0; 0,1,0,0; 0,0,0,0;];
%  F = D_A - D_A'
%  F_indir = TensorEM([0,0,0], [0,0,1])
%  
%  D_A_boost = BoostDerivadas4Vector(D_A, u);
%  F_boost = D_A_boost - D_A_boost'
%  [E_boost, B_boost] = Boost_EM([0,0,0], [0,0,1], u)
%  
%  F_boost_indir = TensorEM(E_boost, B_boost)