function S_uvr = TensorFlujoSpin(Phi, mu, ni, rho) 

matriz = MatrizFlujoSpin(mu, ni, rho);

S_uvr = Phi' * matriz * Phi;

end