% ----------------------    PRUEBAS MAYORES    ----------------------

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

% Estas probaturas se limitan a probar estas cosas:
% ·	Fase
% ·	Módulo global (jt)
% ·	Momento
% ·	Spin origen (en el spinor_A).


clear;

[gt, gx, gy, gz] = MatricesGamma();

% Prueba 1: confirmar que en las soluciones de ondas planas son correctas y
% por tanto el Hamiltoniano al aplicar las derivadas y las matrices
% adecuadas genera un bispinor paralelo al bispinor original (es decir la
% onda con su bispinor y su fase es solución del Hamiltoniano).
% H*Phi = E*Phi

 spinor_base = [1, 0];
%v = [0, 0, 0];             % OK
%v = [0.4, 0, 0];           % OK
%v = [0, 0, 0.4];           % OK
%v = [0, -0.4, 0.55];       % OK
%v = [0.3, -0.5, 0.45];      % OK
%p = fGamma(v) * v;

% E = sqrt(1+norm(p)^2);
% 
% 
% bispinor = DiracSpinorPlainWave(p, spinor_base)
% 
% gv = (p(1)*gx + p(2)*gy + p(3)*gz);
% 
% matriz_mov = gt * gv;
% matriz_t = gt;
% 
% dif_fase = (1i * (matriz_mov + matriz_t) * bispinor')'
% 
% % La fase en el tiempo cambia con la misma forma que la onda.
% dfase_esperado = bispinor * E * (-1i)



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

spinor_base = [1, 0];   %Z

%%% dir_dif_"no_norm", spoiler: siempre estuvo normalizada.
%dir_dif_no_norm = [1,0,0];
%v = [0, 0, 0];
% Resultado: d_vx = -2, lo demás 0.

% dir_dif_no_norm = [0,0,1];
% v = [0, 0, 0];
% Resultado: d_vz = -2, lo demás 0.

% dir_dif_no_norm = [1,0,0];
% v = [0.4, 0, 0];
% Resulado: d_dens_prob =   -0.8000    d_vx = -1.6800

% dir_dif_no_norm = [0,0,1];
% v = [0, 0, 0.4];
% Resulado: d_dens_prob =   -0.8000    d_vz = -1.6800

% dir_dif_no_norm = [1,0,0];
% v = [0, 0, 0.4];
% Resultado: d_vx = -2   d_sx = -0.4174

% dir_dif_no_norm = [0,0,1];
% v = [0.4, 0, 0];
% Resultado: d_vz = -2   d_sx = 0.4174

% dir_dif_no_norm = [1,0,0];
% v = [0, 0.4, 0];
% Resultado: d_vx = -2

% dir_dif = NormalizaCustom(dir_dif_no_norm, 1);
% p = fGamma(v) * v;
% E = sqrt(1+norm(p)^2);
%  
% bispinor = DiracSpinorPlainWave(p, spinor_base);

% La derivada por definición es 1, pero sólo en la dirección especificada
% en la prueba (quiero saber el efecto de haber dicha derivada sobre el
% spin.
% gv = dir_dif(1)*gx + dir_dif(2)*gy + dir_dif(3)*gz;
% 
% dif_amplitud = -(gt * gv * bispinor')';
% 
% [d_dens_prob, d_vx, d_vy, d_vz, d_sx, d_sy, d_sz] = ObtenVariacionPropiedadesIntrinsecas(bispinor, dif_amplitud);
% 
% props_intrinsecas = [d_dens_prob, d_vx, d_vy, d_vz, d_sx, d_sy, d_sz]


% Prueba3 esta vez las variaciones de la onda serán aquellas que hacen que
% la onda en un punto y en otro tengan una diferencia de dp, sea dp_x o
% dp_y o dp_z, esto es más jodido y hace que aparezcan 4 variables en
% juego:
% · la dirección del dp.
% · la dirección de la derivada
% · el spin (que siempre es Z).
% · la dirección de la velocidad (inicialmente 0).

spinor_base = [1, 0];   %Z

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

% dir_dif_no_norm = [0,0,1];
% dp_no_norm = [0, 0, 1];
% v = [0, 0, 0.4];
% Resultado:  d_dens_prob = 0.9165

%dir_dif_no_norm = [0,0,1];
%dp_no_norm = [0, 0, 1];
%v = [0.4, 0, 0];
% Resultado:  d_dens_prob = 1.0000 d_vx = -0.1753

%dir_dif_no_norm = [0,0,1];
%dp_no_norm = [1, 0, 0];
%v = [0, 0, 0.4];
% Resultado: d_vx = -0.1913   d_sx = -0.9564

% dir_dif_no_norm = [0,0,1];
% dp_no_norm = [1, 0, 0];
% v = [0.4, 0, 0];
% Resultado: d_vz = 0.3666    d_sx =  -0.9165

% dir_dif_no_norm = [0,0,1];
% dp_no_norm = [1, 0, 0];
% v = [0, 0.4, 0];
% Resultado: d_sx =   -0.9564

%dir_dif_no_norm = [1,0,0];
%dp_no_norm = [0, 0, 1];
%v = [0, 0, 0.4];
% Resultado:  d_vx = 0.3666      d_sx = 0.9165

%dir_dif_no_norm = [1,0,0];
%dp_no_norm = [1, 0, 0];
%v = [0, 0, 0.4];
% Resultado: d_dens_prob = 0.9165   d_vz = -0.1753

%dir_dif_no_norm = [1,0,0];
%dp_no_norm = [0, 0, 1];
%v = [0.4, 0, 0];
% Resultado: d_dens_prob =  d_vz = -0.1913  d_sx = 0.9564

% dir_dif_no_norm = [1,0,0];
% dp_no_norm = [1, 0, 0];
% v = [0.4, 0, 0];
% Resultado: d_dens_prob =  0.9165

% dir_dif_no_norm = [1,0,0];
% dp_no_norm = [0, 1, 0];
% v = [0.4, 0, 0];
% Resultado: d_vy = -0.1913

% dir_dif_no_norm = [1,0,0];
% dp_no_norm = [0, 1, 0];
% v = [0, 0, 0.4];
% Resultado: todo 0!

% dir_dif_no_norm = [1,0,0];
% dp_no_norm = [0, 1, 0];
% v = [0, 0.4, 0];
% Resultado: d_vx =  0.3666

% dir_dif = NormalizaCustom(dir_dif_no_norm, 1);
% dp = NormalizaCustom(dp_no_norm, 1);
% p = fGamma(v) * v;
% E = sqrt(1+norm(p)^2);
% 
% bispinor = DiracSpinorPlainWave(p, spinor_base);
% 
% dPhi_por_dp = ObtenDiferencialMomento(bispinor, p, dp)

% % Cuando la verificación que viene luego no iba bien, tuve que hacer todo 
% % esto para entender qué estaba haciendo mal... lo mío me costó.
% %bispinor2 = DiracSpinorPlainWave(p + [0,0,0.0001], spinor_base);
% bispinor2 = DiracSpinorPlainWave(p + [0.0001,0,0], spinor_base);
% dif_spinor_unitario = 10000 * (bispinor2 - bispinor)
% dens1 = bispinor * bispinor';           
% dens2 = bispinor2 * bispinor2';
% dif_dens_norm = 10000 * (dens2 - dens1) / dens1;
% factor1 = dPhi_por_dp(1) / dif_spinor_unitario(1);
% factor2 = dPhi_por_dp(3) / dif_spinor_unitario(3);


% A modo de verificación. Hizo falta introducir el factor dos del bra y 
% el ket en ObtenVariacionPropiedadesIntrinsecas para salvar esto!
%[d_dens_prob, d_vx, d_vy, d_vz, d_sx, d_sy, d_sz] = ObtenVariacionPropiedadesIntrinsecas(bispinor, dPhi_por_dp)
% return

%gv = dir_dif(1)*gx + dir_dif(2)*gy + dir_dif(3)*gz;
 
%dPhi = (gt * gv * dPhi_por_dp')';

%[d_dens_prob, d_vx, d_vy, d_vz, d_sx, d_sy, d_sz] = ObtenVariacionPropiedadesIntrinsecas(bispinor, dPhi)



% Prueba 4: introducir rotaciones en el spin.

spinor_base = [1, 0];   %Z

% dir_dif_no_norm = [1,0,0];
% eje_rotacion_no_norm = [0, 0, 1];
% v = [0, 0, 0];
% Resultado: no valido, rotar sobre el eje Z  esto no cambia el spin (que
% ya está en Z, pero sí genera un cambio de fase, por tanto,
% genera un momento y obviamente el momento tiene efecto en
% el hamiltoniano.

%dir_dif_no_norm = [0,0,1];
%eje_rotacion_no_norm = [1, 0, 0];
%v = [0, 0, 0];
% Resultado:  d_vy = -1

% dir_dif_no_norm = [0,1,0];
% eje_rotacion_no_norm = [1, 0, 0];
% v = [0, 0, 0];
% Resultado: d_vz = 1

%dir_dif_no_norm = [1,0,0];
%eje_rotacion_no_norm = [1, 0, 0];
%v = [0, 0, 0];
% Resultado: nada cambia

% dir_dif_no_norm = [0,0,1];
% eje_rotacion_no_norm = [1, 0, 0];
% v = [0, 0, 0.4];
% Resultado: d_vy = -0.9165  d_sy = 0.2087

% dir_dif_no_norm = [1,0,0];
% eje_rotacion_no_norm = [1, 0, 0];
% v = [0, 0, 0.4];
% Resultado: nada cambia si dir_dif_no_norm y eje_rotacion_no_norm paral

% dir_dif_no_norm = [0,1,0];
% eje_rotacion_no_norm = [1, 0, 0];
% v = [0, 0, 0.4];
% Resultado: d_vz = 0.9165

%dir_dif_no_norm = [0,0,1];
%eje_rotacion_no_norm = [1, 0, 0];
%v = [0.4, 0, 0];
% Resultado: d_vy = -1

%dir_dif_no_norm = [0,1,0];
%eje_rotacion_no_norm = [1, 0, 0];
%v = [0.4, 0, 0];
% Resultado: d_vz = 1   d_sx =  -0.2087

%dir_dif_no_norm = [0,0,1];
%eje_rotacion_no_norm = [1, 0, 0];
%v = [0, 0.4, 0];
% Resultado: d_vy = -0.9165

%dir_dif_no_norm = [0,1,0];
%eje_rotacion_no_norm = [1, 0, 0];
%v = [0, 0.4, 0];
% Resultado: d_vz = 0.9165    d_sy = 0.2087


% dir_dif = NormalizaCustom(dir_dif_no_norm, 1);
% eje_rotacion = NormalizaCustom(eje_rotacion_no_norm, 1);
% p = fGamma(v) * v;
% E = sqrt(1+norm(p)^2);
% 
% bispinor = DiracSpinorPlainWave(p, spinor_base);
% 
% dPhi_por_rot = ObtenDiferencialRotacion(bispinor, p, eje_rotacion)

% chequeo de control.
%[d_dens_prob, d_vx, d_vy, d_vz, d_sx, d_sy, d_sz] = ObtenVariacionPropiedadesIntrinsecas(bispinor, dPhi_por_rot)
%return

% gv = dir_dif(1)*gx + dir_dif(2)*gy + dir_dif(3)*gz;
%  
% dPhi = (gt * gv * dPhi_por_rot')';
% 
% [d_dens_prob, d_vx, d_vy, d_vz, d_sx, d_sy, d_sz] = ObtenVariacionPropiedadesIntrinsecas(bispinor, dPhi)



% Prueba 5 - BONUS TRACK: cambios de fase "injustos", es decir que no son
% coherentes con  la onda, por ejemplo, ver qué pasará con una onda cuyo
% spinor tiene corriente 0 pero que tiene cambios de fase, afectará al
% hamiltoniano...
% Añadiré el efecto de la masa para que cuando el cambie de fase sea 0, las
% propiedades intrínsecas no cambien.

spinor_base = [1, 0];   %Z

%%% dir_dif_"no_norm", spoiler: siempre estuvo normalizada.
 %dir_dif_no_norm = [1,0,0];
 %v = [0, 0.4, 0];
% Resultado: d_dv = -2

% dir_dif_no_norm = [0,1,0];
% v = [0, 0.4, 0];
% Resultado: d_dens_prob: -0.8  d_dv = -1.68

% dir_dif_no_norm = [1,0,0];
% v = [0, 0, 0];
% Resultado: d_vy = -2

% dir_dif_no_norm = [0,1,0];
% v = [0, 0, 0];
% Resultado: d_vx = -2

%dir_dif_no_norm = [0,0,1];
%v = [0, 0, 0];
% Resulado: d_vz = -2

% dir_dif = NormalizaCustom(dir_dif_no_norm, 1);
% p = fGamma(v) * v;
% E = sqrt(1+norm(p)^2);
%   
% bispinor = DiracSpinorPlainWave(p, spinor_base);
%  
% gv = dir_dif(1)*gx + dir_dif(2)*gy + dir_dif(3)*gz;
% 
% dif_Phi = -(gt * gv * bispinor')';
% 
% [d_dens_prob, d_vx, d_vy, d_vz, d_sx, d_sy, d_sz] = ObtenVariacionPropiedadesIntrinsecas(bispinor, dif_Phi);
% 
% dif_props_intrinsecas = [d_dens_prob, d_vx, d_vy, d_vz, d_sx, d_sy, d_sz]


% Prueba 6: esta es diferente a la otras: ver que propiedades intrínsecas
% tiene un bispinor al que le he cambiado pi/2 la fase del spinor_B
% Ver también si su energía cambia.

spinor_base = [1, 0];   %Z

% v = [0.4, 0, 0];
% Resultado: la velocidad pasa a ser en dirección -Y en vez de en X.

% v = [0, 0, 0.4];
% Resultado: la velocidad pasa a ser 0.
% Estamos en algo que ya no es un autoestado y por tanto el spinor
% evolucionara, además dio esta energía-momento, una vez definimos
% correctamente la evolución temporal:
% cuatro_momento_trucado = 0.9165         0         0    0.4364
% Esto ha de estar terminantemente prohibido!

%v = [0, 0.4, 0];
v = [0.35, 0.4, -0.5];

p = fGamma(v) * v;
E = sqrt(1+norm(p)^2);
  
bispinor = DiracSpinorPlainWave(p, spinor_base);

%%%bispinor_trucado = [bispinor(1:2), 1i*bispinor(3:4)];

[jt, jx, jy, jz]  = ObtenCorrientesBispinor(bispinor);

corrientes_intrinsecas = [jt, jx, jy, jz]

grad_Phi = ObtenGradienteOndaPlana(bispinor, [E, p]);

[pt, px, py, pz] = ObtenEnergiaMomentoPruebas(bispinor, grad_Phi);

cuatro_momento = [pt, px, py, pz]

%[jt_tr, jx_tr, jy_tr, jz_tr]  = ObtenCorrientesBispinor(bispinor_trucado);
%
% corrientes_trucadas = [jt_tr, jx_tr, jy_tr, jz_tr]
% 
% grad_Phi_trucado = ObtenGradienteOndaPlana(bispinor_trucado, [E, p]);
% 
% [pt_truc, px_truc, py_truc, pz_truc, dPhi_dt_verdadero] = ObtenEnergiaMomentoPruebas(bispinor_trucado, grad_Phi_trucado);
% 
% cuatro_momento_trucado = [pt_truc, px_truc, py_truc, pz_truc]

 
% [d_dens_prob, d_vx, d_vy, d_vz, d_sx, d_sy, d_sz] = ObtenVariacionPropiedadesIntrinsecas(bispinor, dPhi_dt_verdadero);
% 
% variacion_props_intrinsecas_onda_trucada = [d_dens_prob, d_vx, d_vy, d_vz, d_sx, d_sy, d_sz]
 


% dPhi_dt_verdadero es necesario para acabar con el fraude! 
function [pt, px, py, pz] = ObtenEnergiaMomentoPruebas(phi, grad_Phi)
    
    dPhi_dt = grad_Phi(1,:);
    dPhi_dx = grad_Phi(2,:);
    dPhi_dy = grad_Phi(3,:);
    dPhi_dz = grad_Phi(4,:);
    
    dens_prob = norm(phi)^2;
    
    phi_conj = conj(phi);
    
    pt = -1i * (phi_conj * dPhi_dt.')' / dens_prob;
    px = 1i *(phi_conj * dPhi_dx.')' / dens_prob;
    py = 1i *(phi_conj * dPhi_dy.')' / dens_prob;
    pz = 1i *(phi_conj * dPhi_dz.')' / dens_prob;
    
    %Hacer las cosas bien...
    [at, ax, ay, az] = MatricesAlfa();
    gt = MatrizGamma(0); 
    
    % El ultimo termino es el efecto de la masa.
    dPhi_dt_verdadero = ((ax * dPhi_dx.' + ay * dPhi_dy.' + az * dPhi_dz.') + 1i * gt * phi.').';
    % Con este valor está manifiestamente mal, hay que añadir el gt a la masa.
    % dPhi_dt_verdadero = ((ax * dPhi_dx.' + ay * dPhi_dy.' + az * dPhi_dz.') + 1i *  phi.').';     
    
    dif_dt = dPhi_dt_verdadero - dPhi_dt;

    % Debug
%     ratio1 = -i * dPhi_dt_verdadero(1) / phi(1)
%     ratio2 = -i * dPhi_dt_verdadero(2) / phi(2)
%     ratio3 = -i * dPhi_dt_verdadero(3) / phi(3)
%     ratio4 = -i * dPhi_dt_verdadero(4) / phi(4)
    
    pt = -1i * (phi_conj * dPhi_dt_verdadero.') / dens_prob;  
    
end

% Esta es la función que da los resultados importantes en las pruebas 2,3,4
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
    
    [alfa_t, alfa_x, alfa_y, alfa_z] = MatricesAlfa();
    
    dens_prob = Phi' * alfa_t *  Phi;
    
    jx = Phi' * alfa_x * Phi;
    jy = Phi' * alfa_y * Phi;
    jz = Phi' * alfa_z * Phi;
    
    Sx = spinor_A_norm' * PauliX * spinor_A_norm;
    Sy = spinor_A_norm' * PauliY * spinor_A_norm;
    Sz = spinor_A_norm' * PauliZ * spinor_A_norm;
    
    d_dens_prob = (Phi' * alfa_t * dPhi);
    
    d_jx = Phi' * alfa_x * dPhi;
    d_jy = Phi' * alfa_y * dPhi;
    d_jz = Phi' * alfa_z * dPhi;
    
    d_Sx = spinor_A_norm' * PauliX * dPhi_Spinor_A_norm;
    d_Sy = spinor_A_norm' * PauliY * dPhi_Spinor_A_norm;
    d_Sz = spinor_A_norm' * PauliZ * dPhi_Spinor_A_norm;    
    
    d_vx = d_jx / dens_prob - d_dens_prob  * jx / dens_prob^2;
    d_vy = d_jy / dens_prob - d_dens_prob  * jy / dens_prob^2;
    d_vz = d_jz / dens_prob - d_dens_prob  * jz / dens_prob^2;

    d_sx = d_Sx - Sx * dAmplitud_Spinor_A;
    d_sy = d_Sy - Sy * dAmplitud_Spinor_A;
    d_sz = d_Sz - Sz * dAmplitud_Spinor_A;
    
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
    d_amplitud = 1/(2*sqrt(2)) * dE_dp / sqrt(E+1)  / normA;
   
    dif_spinor_A = + d_amplitud * spinor_A';
    
    dif_spinor_B = (dp_dot_sigmas/(E+1) - (p_dot_sigmas * dE_dp) / (E+1)^2) * spinor_A' + d_amplitud * spinor_B';
    
    % Recordar que un cambio en el momento no afecta para nada el spin A.
    dPhi = [dif_spinor_A', dif_spinor_B'];

end

function dPhi = ObtenDiferencialRotacion(bispinor, p, eje_rotacion)
    
    operador_rotacion = 1i * 0.5 * PauliVectorEscalarProd(eje_rotacion);
    
    spinor_A_norm = bispinor(1:2) / norm(bispinor(1:2));
    
    spinor_rotacion = operador_rotacion * spinor_A_norm';
    
    dPhi = DiracSpinorPlainWave(p, spinor_rotacion');

end

function grad_Phi = ObtenGradienteOndaPlana(Phi, p4)

    dPhi_dt = -1i * Phi * p4(1);
    dPhi_dx = 1i * Phi * p4(2);
    dPhi_dy = 1i * Phi * p4(3);
    dPhi_dz = 1i * Phi * p4(4);

    grad_Phi = [dPhi_dt; dPhi_dx; dPhi_dy; dPhi_dz];
    
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
% bispinor = DiracSpinorPlainWave(p, spinor_base);
% 
% [sx_pos, sy_pos, sz_pos] = ObtenSpin(bispinor, 1);
% [sx_neg, sy_neg, sz_neg] = ObtenSpin(bispinor, -1);
% s_pos = [sx_pos, sy_pos, sz_pos]
% s_neg = [sx_neg, sy_neg, sz_neg]
% s = (s_pos + s_neg) / 2
