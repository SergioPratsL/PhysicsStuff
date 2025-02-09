clear;

% OK
% spinor = [1, 0];
% dir_origen = SpinorToVector(spinor)
% eje = [1, 0, 0];
% radianes = pi/4;
% sp_rot = RotacionSpin(spinor, eje, radianes)
% dir_rotada = SpinorToVector(sp_rot')

% spinor = [1, 1];
% dir_origen = SpinorToVector(spinor)
% eje = [1, 0, 0];
% radianes = pi/2;
% sp_rot = RotacionSpin(spinor, eje, radianes)
% dir_rotada = SpinorToVector(sp_rot')

% spinor = [1, 1];
% dir_origen = SpinorToVector(spinor)
% eje = [0 ,1, 1];
% radianes = -pi/2;
% sp_rot = RotacionSpin(spinor, eje, radianes)
% dir_rotada = SpinorToVector(sp_rot')

% spinor = [1, 1];
% dir_origen = SpinorToVector(spinor)
% eje = [0 ,1, 1];
% radianes = pi/3;
% sp_rot = RotacionSpin(spinor, eje, radianes)
% dir_rotada = SpinorToVector(sp_rot')


% Ahora con boosts!

% spinor = [1, 0];
% dir_origen = SpinorToVector(spinor)
% v = [0, 0, 0.4];
% sp_boost = Boost(spinor, v)
% dir_rotada = SpinorToVector(sp_boost')

%spinor = [1, 0];
% dir_origen = SpinorToVector(spinor)
% v = [-0.4, 0, 0];
% sp_boost = BoostSpin(spinor, v)
% dir_rotada = SpinorToVector(sp_boost')


% No me fío de los resultados obtenidos transformando Paulis... a usar
% Diracs...
p0 = [0, 0, 0];

% spin_base = [1, 0];   
% bispinor = DiracSpinorPlainWave(p0, spin_base)
% v = [0, 0, 0.4];
% biSpinor_boost = BoostBispinor(bispinor, v)
% Sigma = fGamma(v);
% biSpinor_esperado = DiracSpinorPlainWave(Sigma*v, spin_base)
% chequeo_espero_sea_0 = biSpinor_boost(3)/biSpinor_boost(1) - biSpinor_esperado(3)/biSpinor_esperado(1)
%Tras los cambios dio un resultado bueno!


% spin_base = [1, 0];   
% bispinor = DiracSpinorPlainWave(p0, spin_base)
% %v = [0.4, 0, 0.4];
% v = [0, 0.4, 0.4];
% bispinor_boost = BoostBispinor(bispinor, v)
% ph_rot = bispinor_boost(1:2);
% dir_principal = SpinorToVector(ph_rot')
% ph_helic = bispinor_boost(3:4);
% dir_helicidad = SpinorToVector(ph_helic')
% Sigma = fGamma(v);
% bispinor_esperado = DiracSpinorPlainWave(Sigma*v, spin_base)
% chequeo_espero_sea_0 = bispinor_boost(4)/bispinor_boost(1) - bispinor_esperado(4)/bispinor_esperado(1);
% dir_helicidad_esperada = SpinorToVector(bispinor_esperado(3:4))
%Si uso -X, Y o -Y como direccion de la velocidad, la helicidad no cambia,
% vx --> bispinor_boost = 1.0518         0    0.2306    0.2306
% vy --> bispinor_boost = 1.0518         0    0.2306   -0.2306i


spin_base = [1, 0];   
bispinor = DiracSpinorPlainWave(p0, spin_base)
%v = [0.3, 0.3, 0.01];
%v = -[0.3, 0.3, 0.01];
v = -[0.6, 0.4, -0.3];
bispinor_boost = BoostBispinor(bispinor, v);
ph_rot = bispinor_boost(1:2);
dir_principal = SpinorToVector(ph_rot')
ph_helic = bispinor_boost(3:4);
dir_helicidad = SpinorToVector(ph_helic')
Sigma = fGamma(v);
bispinor_boost = bispinor_boost.'
bispinor_esperado = DiracSpinorPlainWave(Sigma*v, spin_base)
chequeo_espero_sea_0 = bispinor_boost(4)/bispinor_boost(1) - bispinor_esperado(4)/bispinor_esperado(1);
dir_helicidad_esperada = SpinorToVector(bispinor_esperado(3:4))
% Dos boosts de velocidade opuestas dan lugar a la misma helicidad, pero la
% fase del segundo  spinor es opuesta para cada una de las velocidades.
% con v positiva: 1.0258 + 0.0000i   0.0000 + 0.0000i   0.0054 + 0.0000i   0.1615 + 0.1615i
% con v negativa: 1.0258 + 0.0000i   0.0000 + 0.0000i  -0.0054 - 0.0000i  -0.1615 - 0.1615i

% spin_base = [1, 0];   
% bispinor = DiracSpinorPlainWave(p0, spin_base)
% v = sqrt(2)*[0, 0.2, 0.2];
% biSpinor_boost = BoostBispinor(bispinor, v)
% ph_rot = biSpinor_boost(1:2);
% dir_principal = SpinorToVector(ph_rot')
% ph_helic = biSpinor_boost(3:4);
% dir_helicidad = SpinorToVector(ph_helic')
% Sigma = fGamma(v);
% biSpinor_esperado = DiracSpinorPlainWave(Sigma*v, spin_base)
% chequeo_espero_sea_0 = biSpinor_boost(3)/biSpinor_boost(1) - biSpinor_esperado(3)/biSpinor_esperado(1);
% dir_helicidad_esperada = SpinorToVector(biSpinor_esperado(3:4))
% biSpinor_boost = [1.0225 + 0.0000i   0.0000 + 0.0000i   0.1509 + 0.0000i   0.0000 - 0.1509i]
% La helicidad girará en el plano s-v el doble del ángulo entre p y v...



% spin_base = [1, 0];   
% bispinor = DiracSpinorPlainWave(p0, spin_base)
% eje = [1, 0, 0];
% radianes = pi/2;
% bispinor_rot = RotacionBispinor(bispinor, eje, radianes)
% sp_rot = bispinor_rot(1:2);
% dir_rotada = SpinorToVector(sp_rot')
% Ahora va bien con la fórmula de la wikipedia (no como los boosts).


% p1 = [0, 0.4, 0];
% spin_base = [1, 0];   
% bispinor = DiracSpinorPlainWave(p1, spin_base)
% eje = [1, 0, 0];
% radianes = -pi/2;
% bispinor_rot = RotacionBispinor(bispinor, eje, radianes)
% sp_rot = bispinor_rot(1:2);
% dir_rotada = SpinorToVector(sp_rot')
% ph_helic = bispinor_rot(3:4);
% dir_helicidad = SpinorToVector(ph_helic')
% La rotación no falla...


% Spines muy parecidos para diferentes ondas planas:
% spin_base = [1, 0];   
% v1 = [0.4, 0, 0];
% Sigma1 = fGamma(v1);
% bispinor1 = DiracSpinorPlainWave(Sigma1*v1, spin_base)
% mom4_1 = MomentoEnergiaBispinor(bispinor1)
% veloc_1 = mom4_1(2:4)/mom4_1(1)
% 
% v2 = 0.4*[sqrt(3)/2, .5, 0];
% Sigma2 = fGamma(v2);
% bispinor2 = DiracSpinorPlainWave(Sigma2*v2, spin_base)
% mom4_2 = MomentoEnergiaBispinor(bispinor2)
% veloc_2 = mom4_2(2:4)/mom4_2(1)
% La fase relativa cambia pi/6 entre esto dos bispinors, es decir, 30ª


% v2 = 0.4*[-0.75, .25*sqrt(3), 0.5];
% Sigma2 = fGamma(v2);
% bispinor2 = DiracSpinorPlainWave(Sigma2*v2, spin_base)
% mom4_2 = MomentoEnergiaBispinor(bispinor2)
% veloc_2 = mom4_2(2:4)/mom4_2(1)
% Parece que el operador alfa hace su trabajo muy bien



% Pruebas para ver que la densidad de probabilidad en ondas planas compensa
% la contracción de longitudes.
% spin_base = [1, 0];   
% bispinor_0 = DiracSpinorPlainWave(p0, spin_base)
% 
% v = [0, 0.4, 0];
% bispinor_boost = BoostBispinor(bispinor_0, v)
% 
% gt = MatrizGamma(0);
% 
% dens_v0 = bispinor_0 * gt * gt * bispinor_0'
% dens_v = bispinor_boost' * gt * gt * bispinor_boost
% Fue bien, la densidad de probabilidad crece para compensar la contracción
% de longitudes.



% spin_base = [1, 0];   
% bispinor_boost_neg = DiracSpinorPlainWave(p0, spin_base, -1)
% 
% v = [0.4, 0, 0];
% pv = fGamma(v)*v;
% bispinor_v = DiracSpinorPlainWave(pv, spin_base, -1)
