function Ms = MatrizMs( index1, index2 )
% Incluye el factor -i/4

Gamma1 = MatrizGamma(index1);
Gamma2 = MatrizGamma(index2);

%Ms = -(1i/4) * (Gamma1*Gamma2' - Gamma2*Gamma1');
Ms = -(1i/4) * (Gamma1*Gamma2 - Gamma2*Gamma1);

% ñapa:
%Ms = -(1i/4) * (Gamma1*conj(Gamma2') - Gamma2*conj(Gamma1')); Siempre 0.
%Ms = -(1i/4) * (Gamma1*Gamma2 + Gamma2*Gamma1)
