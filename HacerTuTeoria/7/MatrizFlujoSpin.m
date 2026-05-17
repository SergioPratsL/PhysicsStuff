function matriz = MatrizFlujoSpin(mu, ni, rho) 

gt = MatrizGamma(0);

g_mu = MatrizGamma(mu);
g_ni = MatrizGamma(ni);
g_rho = MatrizGamma(rho);

% Parece que dan los mismos resultados...
%matriz = (1i/2) * g_rho * (g_mu*g_ni - g_ni*g_mu);
matriz = (1i/2) * gt * g_rho * (g_mu*g_ni - g_ni*g_mu);

end