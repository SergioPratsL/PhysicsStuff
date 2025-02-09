% Este es el ejemplo más trivial de obtener el mismo resultado que la
% transformación de velocidades con muchas iteraciones, vBC_actual es una
% variable de control, lo que importa es el resultado de  vAC_actual.
% La razón de ello es que los SRIs de A y C están alineados con B por lo
% que no hay problemas entre ellos, pero sí entre A y C (que deben rotar).
% En este ejemplo se calcula cómo cambia la velocidad de C vista por A
% haciendo pasar la aceleración propia (dif_vBC) primero a B, medienta
% boost de -vBC_actual y luego a C con un boost de vBC.
% Es el ejemplo trivial
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
pAC_actual = vAC_actual*SigmaAC_actual;

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
    
    %dif_pAC_actual = Boost(dif_p_local, -vAC_actual);
    dif_pAC_actual = Boost(dif_pBC_actual, vBA);
    pAC_actual = pAC_actual + dif_pAC_actual(2:4);
    SigmaAC_actual = sqrt(1+norm(pAC_actual)^2);
    vAC_actual = pAC_actual / SigmaAC_actual; 
    
    n = n + 1;
end

n_fin = n
vBC = vBC_actual
vAC = vAC_actual

vAC_transf = Vel_Addition_Law(vBC_final, vBA)

% Fue bien :)
