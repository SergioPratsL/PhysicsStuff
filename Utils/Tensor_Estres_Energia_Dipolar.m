% �Ojo! Este tensor de estres energia lo he derivado yo, bas�ndome en la
% suposici�n que dM^(mu, alfa)/dx_alfa = 0 (para cualquier valor de mu).
% Esto est� respaldado en las ecuaciones de Maxwell, que l�gicamente no
% aplican al tensor de momentos dipolares M, en el caso de M creo que
% equivale a pedir que las bound currents sean 0.
% Esto es controvertido...
function T = Tensor_Estres_Energia_Dipolar(F, P)
  eta = MatrizMinkowski();
  
  % Caso original del EM que va bien.
  %T = -(F*eta*F' - 1/4 * eta*trace(F*eta*F'*eta)); 
  
  % Este es el valor que da la energ�a correcta.
  T =  (F*eta*P' + P*eta*F') - 1/4 * eta * (trace(F*eta*P'*eta) + trace(P*eta*F'*eta)); 
end