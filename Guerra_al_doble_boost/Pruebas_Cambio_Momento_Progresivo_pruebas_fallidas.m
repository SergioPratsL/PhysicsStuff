clear;

% Prueba 1
vBA = [0.4, 0, 0];
vBC_final = [0, 0.4, 0];
% Harán falta algo más de 4000 iteraciones (4255)
dif_vBC = [0, 0.0001, 0];
dif_p_local = [0, dif_vBC];
% % Resultado: vAC = [-0.3853    0.3688         0];

%Prueba 1 parte 2, ahora A hace de B y B hace de C, con lo que C hace de A
% vBA = [-0.3853,0.3688,  0];
% vBC_final = [-0.4, 0, 0];
% dif_vBC = [-0.0001, 0, 0];
% dif_p_local = [0, dif_vBC];
% vAC = [0.0165   -0.3842         0]


% Esto siempre tiene que ser así!
vBC_actual = [0,0,0];
vAC_actual = -vBA;

SigmaBC_actual = fGamma(vBC_actual);
SigmaAC_actual = fGamma(vAC_actual);
pBC_actual = vBC_actual*SigmaBC_actual;
pAC_actual = SigmaAC_actual * [1, vAC_actual];

check_ini = pAC_actual(1)^2 - norm(pAC_actual(2:4))^2

% Hacemos un boost directo de la diferencia de momento y vamos actualizando
% la velocidad a cada iteración.

vAC_old = [0,0,0];

n = 1;
% La condición del bucle tiene que cambiar de prueba en prueba!
while vBC_actual(2) < vBC_final(2) && n < 10000
%while vBC_actual(1) > vBC_final(1) && n < 10000
    dif_pBC_actual = Boost(dif_p_local, -vBC_actual);
    pBC_actual = pBC_actual + dif_pBC_actual(2:4);
    SigmaBC_actual = sqrt(1+norm(pBC_actual)^2);
    vBC_actual = pBC_actual / SigmaBC_actual;

    if n == 2547
        n = n + 1;
    end
    
    vBA_AddLaw = Vel_Addition_Law( vAC_actual, vBC_actual );
    vAB_AddLaw = Vel_Addition_Law( vBC_actual, vAC_actual );
    eje_rot = cross(vAC_actual, vBC_actual);
    eje_rot = eje_rot / norm(eje_rot);    
    val = -dot(vBA_AddLaw, vAB_AddLaw) / norm(vBA_AddLaw)^2;
    if( abs(val) > 1)
        val = sign(val);
    end
    angulo = acos(val);
    rotMatrix = RotationMatrixGeneral(eje_rot, angulo);
    
    vAC_actual_rot = vAC_actual * rotMatrix';
    dif_p_local_rot = [dif_p_local(1), dif_p_local(2:4) * rotMatrix];
    %dif_pAC_actual = Boost(dif_p_local, -vAC_actual_rot);
    dif_pAC_actual = Boost(dif_p_local_rot, -vAC_actual_rot);
    %dif_pAC_actual = Boost(dif_p_local_rot, -vAC_actual);
    
    %dif_pAC_actual = Boost(dif_p_local, -vAC_actual);
    
    % Coge los valores de la iteración anterior, suficientemente bueno
    %dv_AC = vAC_actual - vAC_old;
    %%%dv_AC = dif_pAC_actual(2:4);     % mal
    
    % Esto falla y casi no hace nada
    %dif_pAC_actual_rot3 = ThomasPrecession(SigmaAC_actual, dv_AC, vAC_actual, dif_pAC_actual(2:4));
    
    % Compara iteraciones diferentes pero el impacto es poco.
%     vBA_AddLaw = Vel_Addition_Law( vAC_actual, vBC_actual );
%     vAB_AddLaw = Vel_Addition_Law( vBC_actual, vAC_actual );
%     eje_rot = cross(vAC_actual, vBC_actual);
%     eje_rot = eje_rot / norm(eje_rot);    
%     val = -dot(vBA_AddLaw, vAB_AddLaw) / norm(vBA_AddLaw)^2;
%     if( val > 1)
%         val = 1;
%     end
%     if( val < -1 )
%         val = -1;
%     end    
%     angulo = acos(val);
%     rotMatrix = RotationMatrixGeneral(eje_rot, angulo);
%     dif_pAC_actual(2:4) = dif_pAC_actual(2:4) * rotMatrix;
    
    pAC_actual = pAC_actual + dif_pAC_actual;
    % Esto también falla.
    %pAC_actual =  ThomasPrecession(SigmaAC_actual, dv_AC, vAC_actual, pAC_actual);
    
    SigmaAC_actual = sqrt(1+norm(pAC_actual(2:4))^2);
    
    %vAC_old = vAC_actual;    
    vAC_actual = pAC_actual(2:4) / SigmaAC_actual; 
    
    n = n + 1;
end

n_fin = n
vBC = vBC_actual;
vAC = vAC_actual

% Usemos los valores del algoritmo (nada cambia)
%vAC_AddLaw = Vel_Addition_Law(vBC_final, vBA)
vAC_AddLaw = Vel_Addition_Law(vBC, vBA);

Check = pAC_actual(1)^2 - norm(pAC_actual(2:4))^2

% Esto va mal, muy mal...



