function Pot_F = LarmorCovarianteMomentum( v, a_prop, dt )
% Al tener la aceleracion propia esto es mejor que usar la fórmula de
% FitzPatrick.
% La energía de origen se multiplica por Sigma, pero el tiempo de origen se
% es ds = dt/Sigma, asi que en verdad nada cambia :).
   
Pot = (2/3) * (a_prop(1)^2+a_prop(2)^2+a_prop(3)^2);
    
% Por ser un boost de un vector original de tipo energia [1,0,0,0].
% [1,0,0,0] --> Sigma * [1, v]    
Pot_F = Pot * dt * [1, v];
    
end

