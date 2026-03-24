function matriz = MatrizFlujoSpin(mu, ni, rho) 

g_mu = MatrizGamma(mu);
g_ni = MatrizGamma(ni);
g_rho = MatrizGamma(rho);

matriz = (1i/2) * g_rho * (g_mu*g_ni - g_ni*g_mu);

end