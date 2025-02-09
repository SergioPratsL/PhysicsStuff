% En este ejemplo verifico que se puede obtener la velocidad final vAC
% esperada a partir de ir acelerando C y evaluando en cada iteración la
% velocidad de A, siendo que A y C están "atados" al SRI de B.

% Esta es la manera más complicada de hacerlo.
% Por un lado hace falta rotar dif_vBC, 
%la variación de velocidad que es medida para B (debido a 
% que es imposible que C se alinee con los ejes de A y b a la vez.

% Por otro lado, dado que C está acelerando (al cambiar de SRI), también
% es necesario aplicar la rotación de Thomas con lo que la nueva velocidad
% obtenida en cada iteración debe rotarse un poco...

clear;

vBA = [0.4, 0, 0];
% Sólo se usa la segunda componente de vBC_final para comparar!
vBC_final = [0, 0.4, 0];
% Harán falta algo más de 4000 iteraciones (4255)
dif_vBC = [0, 0.0001, 0];
dif_p_local = [0, dif_vBC];


% Esto siempre tiene que ser así!
vBC_actual = [0,0,0];
vAC_actual = -vBA;
vAC_old = vAC_actual;
%vAC_old = vAC_actual;

SigmaBC_actual = fGamma(vBC_actual);
SigmaAC_actual = fGamma(vAC_actual);
pBC_actual = vBC_actual*SigmaBC_actual;
pAC_actual = SigmaAC_actual * [1, vAC_actual];

check_ini = pAC_actual(1)^2 - norm(pAC_actual(2:4))^2

dv_lab = [0,0,0];

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

    % Calcular la velocidad haciendo un mini boost a vCnewA!
    % Cómo, si C tuviera el mando.
    
    if n == 2000
        n = 2000;
    end
    
    rotMatrix = GetThomasRotMatrix(vBA, vBC_actual);
    dif_vBC_rot = dif_vBC * rotMatrix;
    
    vect_old_CA = SigmaAC_actual * [1, - vAC_actual];
    %vect_boost_CA = Boost(vect_old_CA, dif_vBC);
    vect_boost_CA = Boost(vect_old_CA, dif_vBC_rot);
    SigmaAC_actual = vect_boost_CA(1);
    
    vAC_actual = - vect_boost_CA(2:4) / SigmaAC_actual;

    % Diferencia de velocidades no local sino vista desde A que hace de lab
    dv_lab = vAC_actual - vAC_old;
    
    vAC_actual_rot = ThomasPrecession(SigmaAC_actual, dv_lab, vAC_old, vAC_actual);
    
    vAC_old =  vAC_actual; % vAC_actual todavia no ha cambiado
    vAC_actual = vAC_actual_rot;
    
    n = n + 1;
end

n_fin = n;
vBC = vBC_actual;
vAC = vAC_actual

% Usemos los valores del algoritmo (nada cambia)
%vAC_AddLaw = Vel_Addition_Law(vBC_final, vBA)
vAC_AddLaw = Vel_Addition_Law(vBC, vBA)

Check = pAC_actual(1)^2 - norm(pAC_actual(2:4))^2;

prod_esc = dot(vAC_AddLaw, vAC) / norm(vAC_AddLaw) / norm(vAC)


