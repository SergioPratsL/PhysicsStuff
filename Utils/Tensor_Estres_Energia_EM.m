% Contribución decisiva de GeminiIA :D
function T = Tensor_Estres_Energia_EM(F)
  eta = MatrizMinkowski();
    
  T = -(F*eta*F' - 1/4 * eta*trace(F*eta*F'*eta)); 
end