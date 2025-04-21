function SigmaMayusc = MatrizSpin_4_4(index)
% index = 1 para X, 2 para Y, 3 para Z.

pauli = MatrizPauli(index);
m_cero = [0,0;0,0];

SigmaMayusc = [pauli, m_cero; m_cero, pauli];







