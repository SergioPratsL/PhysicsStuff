% Estas pruebas servirán para tener una idea de cómo puede afectar en el
% hamiltoniano el hecho que hayan variaciones en módulo o en la orientación
% del spin o el momento... el bispinor tiene muchos grados de libertad y
% quiero saber cómo afecta su derivada (en la direccion X, por elegir una)
% al spin de un punto.


% Lo primero es verificar que la onda plana genera autoestados, parece
% obvio, pero es mejor empezar por lo más sencillo.

% Sea la exponencial de la onda plana -i(p_u * x_u), si la velocidad es v*e1  
% (siendo e1 la direcció y v un número), entonces tendremos que
% d_e1(Phi) = -i*p*Phi y las derivadas respecto a las direcciones 
% espaciales perpendiculares serían 0.
% A esta derivada hay que aplicarle el correspondiente alfa_e1

% Luego está la derivada temporal, que en ausencia de potencial es
% simplemente m^c2, es decir, 1 y luego estaría el potencial V al que
% tampoco se le aplica ninguna Gamma_t ni nada, porque si no romperia los
% autovectores en las ondas planas.


% Así pues lo primero es ver cómo evoluciona temporalmente la onda plana en
% un punto, luego lo extenderé a cambios de amplitud en las ondas planas a
% ver qué pasa.

clear;

[gt, gx, gy, gz] = MatricesGamma();

% Prueba 1: confirmar que en las soluciones de ondas planas son correctas y
% por tanto el Hamiltoniano al aplicar las derivadas y las matrices
% adecuadas genera un bispinor paralelo al bispinor original (es decir la
% onda con su bispinor y su fase es solución del Hamiltoniano).
% H*Phi = E*Phi

% spin_base = [1, 0];   
% %v = [0, 0, 0];             % OK
% %v = [0.4, 0, 0];           % OK
% %v = [0, 0, 0.4];           % OK
% %v = [0, -0.4, 0.55];       % OK
% v = [0.3, -0.5, 0.45];      % OK
% p = fGamma(v) * v;
% 
% E = sqrt(1+norm(p)^2);
% 
% 
% bispinor = DiracSpinorPlainWave(p, spin_base)
% 
% gv = p(1)*gx + p(2)*gy + p(3)*gz;
% 
% matriz_mov = gt * gv;
% matriz_t = gt;
% 
% dif_fase = (1i * (matriz_mov + matriz_t) * bispinor')'
% 
% % La fase en el tiempo cambia con la misma forma que la onda.
% dfase_esperado = bispinor * E * (-1i)
% 
% % Para comproba que todos los valores dan cero.
% [d_dens_prob, d_vx, d_vy, d_vz, d_sx, d_sy, d_sz] = ObtenVariacionPropiedadesIntrinsecas(bispinor, dif_fase)
% % Con v = [0.3, -0.5, 0.45]; todos los valores daban cero o negligible o
% % imaginario.



% Prueba2: examinar el efecto en el spinor del una diferencia de amplitud
% de fase el cierta dirección, es decir, dPhi_dx = cte * Phi, siendo la
% constante real, por tanto no un cambio de fase, el término gt*m es
% enteramente imaginario por lo que contribuye a la fase pero aquí no hará
% nada. Obtendremos variaciones que no son paralelas a la onda, entonces he
% de ver cómo estas variaciones afectan a las propiedades del spin, es
% decir a la corriente y a la orientación del spin, y al módulo también,
% claro.

% En los resultados ignoro la componente imaginaria de la variación porque
% asumo que al aplicarse sobre el ket se compensaría con la parte que se
% asgina al bra.

spin_base = [1, 0];   %Z

%%% dir_dif_"no_norm", spoiler: siempre estuvo normalizada.
% dir_dif_no_norm = [1,0,0];
% v = [0, 0, 0];
% Resultado: d_vx = 2, lo demás 0.

% dir_dif_no_norm = [0,0,1];
% v = [0, 0, 0];
% Resultado: d_vz = 2, lo demás 0.

% dir_dif_no_norm = [1,0,0];
% v = [0.4, 0, 0];
% Resulado: d_vx = 2.1822

% dir_dif_no_norm = [0,0,1];
% v = [0, 0, 0.4];
% Resultado: d_vz = 2.1822

% dir_dif_no_norm = [1,0,0];
% v = [0, 0, 0.4];
% Resultado: d_vx = 2.1822   d_sx = 0.4174

% dir_dif_no_norm = [0,0,1];
% v = [0.4, 0, 0];
% Resultado: d_vz = 2.1822   d_sx = -0.4174

% dir_dif_no_norm = [1,0,0];
% v = [0, 0.4, 0];
% Resultado: d_vx = 2.1822

% dir_dif = NormalizaCustom(dir_dif_no_norm, 1);
% p = fGamma(v) * v;
% E = sqrt(1+norm(p)^2);
%  
% bispinor = DiracSpinorPlainWave(p, spin_base);

% La derivada por definición es 1, pero sólo en la dirección especificada
% en la prueba (quiero saber el efecto de haber dicha derivada sobre el
% spin.
% gv = dir_dif(1)*gx + dir_dif(2)*gy + dir_dif(3)*gz;
% 
% dif_amplitud = (gt * gv * bispinor')';
% 
% [d_dens_prob, d_vx, d_vy, d_vz, d_sx, d_sy, d_sz] = ObtenVariacionPropiedadesIntrinsecas(bispinor, dif_amplitud)



% Prueba3 esta vez las variaciones de la onda serán aquellas que hacen que
% la onda en un punto y en otro tengan una diferencia de dp, sea dp_x o
% dp_y o dp_z, esto es más jodido y hace que aparezcan 4 variables en
% juego:
% · la dirección del dp.
% · la dirección de la derivada
% · el spin (que siempre es Z).
% · la dirección de la velocidad (inicialmente 0).

spin_base = [1, 0];   %Z

% dir_dif_no_norm = [0,0,1];
% dp_no_norm = [0, 0, 1];
% v = [0, 0, 0];
% Resultado: d_dens_prob = 1

% dir_dif_no_norm = [1,0,0];
% dp_no_norm = [0, 0, 1];
% v = [0, 0, 0];
% Resultado: d_sx = 1

% dir_dif_no_norm = [0,0,1];
% dp_no_norm = [1, 0, 0];
% v = [0, 0, 0];
% % Resultado: d_sx = -1

% dir_dif_no_norm = [1,0,0];
% dp_no_norm = [1, 0, 0];
% v = [0, 0, 0];
% Resultado: d_dens_prob = 1

%dir_dif_no_norm = [1,0,0];
%dp_no_norm = [0, 1, 0];
%v = [0, 0, 0];
% Resultado: Nada cambia

%dir_dif_no_norm = [0,0,1];
%dp_no_norm = [0, 0, 1];
%v = [0, 0, 0.4];
% Resultado:  d_dens_prob = 0.9165

%dir_dif_no_norm = [0,0,1];
%dp_no_norm = [0, 0, 1];
%v = [0.4, 0, 0];
% Resultado:  d_dens_prob = 1.0000 d_vx = -0.2277

%dir_dif_no_norm = [0,0,1];
%dp_no_norm = [1, 0, 0];
%v = [0, 0, 0.4];
% Resultado: d_vx = -0.2087   d_sx = -0.9564

dir_dif_no_norm = [0,0,1];
dp_no_norm = [1, 0, 0];
v = [0.4, 0, 0];
% Resultado: 

dir_dif = NormalizaCustom(dir_dif_no_norm, 1);
dp = NormalizaCustom(dp_no_norm, 1);
p = fGamma(v) * v;
E = sqrt(1+norm(p)^2);

bispinor = DiracSpinorPlainWave(p, spin_base);

dPhi_por_dp = ObtenDiferencialMomento(bispinor, p, dp)

% % Cuando la verificación que viene luego no iba bien, tuve que hacer todo 
% % esto para entender qué estaba haciendo mal... lo mío me costó.
% %bispinor2 = DiracSpinorPlainWave(p + [0,0,0.0001], spin_base);
bispinor2 = DiracSpinorPlainWave(p + [0.0001,0,0], spin_base);
dif_spinor_unitario = 10000 * (bispinor2 - bispinor)
dens1 = bispinor * gt *  bispinor';            % Ouch!!!!
dens2 = bispinor2 * gt * bispinor2';
dif_dens_norm = 10000 * (dens2 - dens1) / dens1;
factor1 = dPhi_por_dp(1) / dif_spinor_unitario(1)
factor2 = dPhi_por_dp(3) / dif_spinor_unitario(3)
 

% A modo de verificación. Hizo falta introducir el factor dos del bra y 
% el ket en ObtenVariacionPropiedadesIntrinsecas para salvar esto!
[d_dens_prob, d_vx, d_vy, d_vz, d_sx, d_sy, d_sz] = ObtenVariacionPropiedadesIntrinsecas(bispinor, dPhi_por_dp)
 return

gv = dir_dif(1)*gx + dir_dif(2)*gy + dir_dif(3)*gz;
 
dPhi = (gt * gv * dPhi_por_dp')';

[d_dens_prob, d_vx, d_vy, d_vz, d_sx, d_sy, d_sz] = ObtenVariacionPropiedadesIntrinsecas(bispinor, dPhi)





function [d_dens_prob_unitaria, d_vx, d_vy, d_vz, d_sx, d_sy, d_sz] = ObtenVariacionPropiedadesIntrinsecas(Phi_ori, dPhi_ori)
    [gt, gx, gy, gz] = MatricesGamma();
   
    Phi = Phi_ori';
    dPhi = dPhi_ori';
    
    % [Sigma_x, Sigma_y, Sigma_z]  = OperadoresSpin_bispinor();
    [PauliX, PauliY, PauliZ] = MatricesPauli();
    spinor_A = Phi(1:2);
    spinor_A_norm = spinor_A / norm(spinor_A);
    dPhi_Spinor_A = dPhi(1:2);
    % Ojo! Se ha de normalizar respecto a la onda, no a la derivada, porque
    % sería aplicar la derivada de la parte proporcional de la onda... este
    % factor de normalización es la muerte.
    dPhi_Spinor_A_norm = dPhi_Spinor_A / norm(spinor_A);
    dAmplitud_Spinor_A = spinor_A_norm' * dPhi_Spinor_A_norm;
    
    alfa_x = gt * gx;
    alfa_y = gt * gy;
    alfa_z = gt * gz;
    
    dens_prob = Phi' * gt *  Phi;
    
    jx = Phi' * alfa_x * Phi;
    jy = Phi' * alfa_y * Phi;
    jz = Phi' * alfa_z * Phi;
    
    Sx = spinor_A_norm' * PauliX * spinor_A_norm;
    Sy = spinor_A_norm' * PauliY * spinor_A_norm;
    Sz = spinor_A_norm' * PauliZ * spinor_A_norm;
%     Sx = Phi' * Sigma_x * Phi;
%     Sy = Phi' * Sigma_y * Phi;
%     Sz = Phi' * Sigma_z * Phi;
    
    d_dens_prob = (Phi' * gt * dPhi);
    
    d_jx = Phi' * alfa_x * dPhi;
    d_jy = Phi' * alfa_y * dPhi;
    d_jz = Phi' * alfa_z * dPhi;
    
    d_Sx = spinor_A_norm' * PauliX * dPhi_Spinor_A_norm;
    d_Sy = spinor_A_norm' * PauliY * dPhi_Spinor_A_norm;
    d_Sz = spinor_A_norm' * PauliZ * dPhi_Spinor_A_norm;    
%     d_Sx = Phi' * Sigma_x * dPhi;
%     d_Sy = Phi' * Sigma_y * dPhi;
%     d_Sz = Phi' * Sigma_z * dPhi;
    
    d_vx = d_jx / dens_prob - d_dens_prob  * jx / dens_prob^2;
    d_vy = d_jy / dens_prob - d_dens_prob  * jy / dens_prob^2;
    d_vz = d_jz / dens_prob - d_dens_prob  * jz / dens_prob^2;

    d_sx = d_Sx - Sx * dAmplitud_Spinor_A;
    d_sy = d_Sy - Sy * dAmplitud_Spinor_A;
    d_sz = d_Sz - Sz * dAmplitud_Spinor_A;
%     d_sx = d_Sx / dens_prob - d_dens_prob * Sx / dens_prob^2;
%     d_sy = d_Sy / dens_prob - d_dens_prob * Sy / dens_prob^2;
%     d_sz = d_Sz / dens_prob - d_dens_prob * Sz / dens_prob^2;
    
    % Es la variación de densidad de probabilidad por unidad de densidad
    % de probabilidad, es decir, en proporción a cuanta densidad de
    % probabilidad tengamos.    
    d_dens_prob_unitaria = d_dens_prob / dens_prob;
    
    
    % Creo que todo se ha de multiplicar por dos y aplicar sólo la parte
    % real, ya que el cambio afectará igual al bra y al ket.
    d_vx = 2 * real(d_vx);
    d_vy = 2 * real(d_vy);
    d_vz = 2 * real(d_vz);
    
    d_sx = 2 * real(d_sx);
    d_sy = 2 * real(d_sy);
    d_sz = 2 * real(d_sz);    
    
    d_dens_prob_unitaria = 2 * real(d_dens_prob_unitaria);
end


function dPhi = ObtenDiferencialMomento(bispinor, p, dp)

    spinor_A = bispinor(1:2);
    spinor_B = bispinor(3:4);
    
    E = sqrt(1+norm(p)^2);
    
    normA = norm(spinor_A);
    
    dE_dp = dot(p, dp) / E;
    
    dp_dot_sigmas = PauliVectorEscalarProd(dp);
    
    p_dot_sigmas = PauliVectorEscalarProd(p);
    
    
    % Debido a que la contracción aumenta, este factor se ha de neutralizar
    % para que el resultado total no cambie dens prob-
    d_amplitud = 1/sqrt(2) * dE_dp / sqrt(E+1) / 2 / normA;
   
    dif_spinor_A = + d_amplitud * spinor_A';
    
    dif_spinor_B = (dp_dot_sigmas/(E+1) - (p_dot_sigmas * dE_dp) / (E+1)^2) * spinor_A' + d_amplitud * spinor_B';
    
    % Recordar que un cambio en el momento no afecta para nada el spin A.
    dPhi = [dif_spinor_A', dif_spinor_B'];

end


% % Prueba para evaluar el spin (es una prueba para saber cómo hacer este script de probaturas).
% Conclusion: ninguna opción es inalterable ante el boost.
% function [sx, sy, sz] = ObtenSpin(Phi, negativo)
%     [Sigma_x, Sigma_y, Sigma_z]  = OperadoresSpin_bispinor(negativo);
%     
%     sx = Phi * Sigma_x * Phi';
%     sy = Phi * Sigma_y * Phi';
%     sz = Phi * Sigma_z * Phi';
% end
% 
% bispinor = DiracSpinorPlainWave(p, spin_base);
% 
% [sx_pos, sy_pos, sz_pos] = ObtenSpin(bispinor, 1);
% [sx_neg, sy_neg, sz_neg] = ObtenSpin(bispinor, -1);
% s_pos = [sx_pos, sy_pos, sz_pos]
% s_neg = [sx_neg, sy_neg, sz_neg]
% s = (s_pos + s_neg) / 2
