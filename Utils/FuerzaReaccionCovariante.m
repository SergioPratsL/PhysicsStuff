function Fr = FuerzaReaccionCovariante(dt, v, a, da_dt, Sigma)
% dt es el intervalo de tiempo medido desde el laboratorio
% v, a y da_dt son magnitudes medidas todas desde el laboratorio

if( Sigma == 0)
    Sigma = fGamma(v);
end

Sigma2 = Sigma^2;
Sigma4 = Sigma2^2;
Sigma6 = Sigma2^3;

dot_v_a = dot(v,a);

term1 = Sigma2 * da_dt;

term2 = Sigma4 * v * dot(v, da_dt);

term3 = 3 * Sigma4 * a * dot_v_a;

term4 = 3 * Sigma6 * v * dot_v_a^2;

Fr = 2/3 * (term1 + term2 + term3 + term4);

