function FR = FR_Covariante( v, a, da, Sigma )
% Esta es la función que calcula la fuerza de reacción covariante... si le
% pasas la velocidad, aceleración y derivada temporal de la aceleración,
% que como no las puedes tener de ya tienes que coger las de la iteración
% anterior... es lo mejor que puedo hacer.
% La velocidad de la luz debe estar normalizada.

    if( Sigma == 0)
        Sigma = fGamma(norm(v));
    end
    
    Sigma2 = Sigma^2;
    Sigma4 = Sigma2^2;
    Sigma5 = Sigma4 * Sigma;

    V_u = Sigma * [1, v];
    
    prod_v_a = dot(v,a);
    prod_v_da = dot(v, da);
    a_cuad = dot(a,a);
    
    A_u_t = Sigma4 * prod_v_a;
    A_u_x = A_u_t * v + Sigma2 * a;
    A_u = [A_u_t, A_u_x];
    
    dtA_u_t = Sigma5 * (4*Sigma2*prod_v_a^2 + prod_v_da + a_cuad);
    dtA_u_x = dtA_u_t * v;
    dtA_u_x = dtA_u_x + (Sigma5*prod_v_a) * a; 
    dtA_u_x = dtA_u_x + 2 * Sigma5 * prod_v_a * a + (Sigma2*Sigma) * da;
    
    % Por fin ensamblamos el jodidísimo 4 vector.
    dA_u = [dtA_u_t, dtA_u_x];

    %FR = (2/3) * (dA_u - (A_u(1)^2 - A_u(2)^2 - A_u(3)^2 - A_u(4)^2) * V_u);
    FR = (2/3) * (dA_u - (-A_u(1)^2 +A_u(2)^2 + A_u(3)^2 + A_u(4)^2) * V_u);

end

