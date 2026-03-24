function DSp = DiracSpinorPlainWave_GtGt(p, Sp2)
% Esto da las ondas planas asumiendo que a la masa también le pertoca un
% término gt de forma que tiene la misma matriz que la energía y el
% potencial EM.

% https://en.wikipedia.org/wiki/Dirac_spinor

% p sería el momento de la onda y phi un two-spinor que apuntará en cierta
% dirección.

if (norm(p) == 0)
    DSp = [[0,0], Sp2];
    return
end

% Primero hay que recuperar Phi de Sp2.

p_square = norm(p)^2;

P_Dot_PauliVector = PauliVectorEscalarProd(p);

Phi_ori = PauliVectorEscalarProd(p) * Sp2.' / p_square;

Sp2_control = PauliVectorEscalarProd(p) * Phi_ori;

E = sqrt(norm(p)^2 + 1);

%%p_dot_fi = P_Dot_PauliVector * fi.';

phi = (E-1)*Phi_ori;
chi = Sp2_control;

DSp = [phi; chi].';

Control_1 = (E-1)*phi - P_Dot_PauliVector*chi
Control_2 = (E-1)*chi - P_Dot_PauliVector*phi

DSp = DSp / norm(DSp);



