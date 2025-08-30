function E = TensorEstresEnergia_Inercial(densidad_masa, v)
% Esto es el tensor de estrés energía de una partícula que se mueve sin
% ejercer ningún traspaso de momento más allá que el momento viaja con la
% partícula

% Densidad masa es muy fullero, no es una masa en reposo, es más bien la
% densidad de probabilidad multiplicada por la masa en reposo, esto es
% importante porque la partícula se contraerá si no está en su
% SRI original.

gamma = fGamma(v);
p = gamma * v;

E = zeros(4);

E(1,1) = gamma;

E(1,2) = p(1);
E(2,1) = p(1);

E(1,3) = p(2);
E(3,1) = p(2);

E(1,4) = p(3);
E(4,1) = p(3);

E(2,2) = p(1)*v(1);
E(3,3) = p(2)*v(2);
E(4,4) = p(3)*v(3);

E(2,3) = p(1)*v(2);
E(3,2) = p(1)*v(2);

E(2,4) = p(1)*v(3);
E(4,2) = p(1)*v(3);

E(3,4) = p(2)*v(3);
E(4,3) = p(2)*v(3);

E = E * densidad_masa;

end