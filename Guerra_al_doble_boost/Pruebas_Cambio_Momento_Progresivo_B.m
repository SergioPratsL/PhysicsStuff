% En este ejemplo se intenta obtener el mismo resultado que la
% transformación de velocidades mediante un algoritmo iterativo evitando
% que la aceleración propia de C tenga que pasar por B, sino enviándola
% directamente a A, que tiene un control de la velocidad de C en todo
% momento, aceptando eso sí, que A está alineado con B por lo que debe
% mediar una rotación entre A y C. Debido a esta rotación, los cambios de
% velocidad propios de C deben rotarse para estar alineados con A y luego
% boostearse con -vAC_actual.
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
    
    rotMatrix = GetThomasRotMatrix(vBA, vBC_actual);
    
    dif_p_local_rot = [dif_p_local(1), dif_p_local(2:4) * rotMatrix];
    dif_pAC_actual = Boost(dif_p_local_rot, -vAC_actual);
    
    pAC_actual = pAC_actual + dif_pAC_actual;
    
    SigmaAC_actual = sqrt(1+norm(pAC_actual(2:4))^2);    
    vAC_actual = pAC_actual(2:4) / SigmaAC_actual; 
    
    n = n + 1;
end

n_fin = n;
vBC = vBC_actual;
vAC = vAC_actual

% Usemos los valores del algoritmo (nada cambia)
%vAC_AddLaw = Vel_Addition_Law(vBC_final, vBA)
vAC_AddLaw = Vel_Addition_Law(vBC, vBA);

Check = pAC_actual(1)^2 - norm(pAC_actual(2:4))^2
