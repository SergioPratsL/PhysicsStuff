clear;

% Espero poder mantener estables estas 2 constantes:
v_ini = [0,0,0]

% Esto es aceleración en partes de la velocidad de la luz por cada unidad
% de tiempo entera, cada iteracion es una milésima de unidad de tiempo.
a_prop_ini = [0, 0, 0];
a_lab_ini = a_prop_ini;     % Verdadero al ser v_ini = 0

% Inicializamos los valores acumulativos
uP_reac = [0,0,0,0];
uP_Larmor = [0,0,0,0];

v = v_ini;
a_prop = a_prop_ini;
a_lab = a_lab_ini;


PR_ini = PotencialReaccion( v_ini, a_prop_ini );

% PRUEBA 1: Aceleración crece y decrece
%% Fase en la que la aceleración crece
% dt_iter = 2.5*10^-4;
% num_iter = 2000;
% daprop_diter = [0.8, 0, 0]/num_iter;
% 
% [v, a_lab, uP_reac_tramo, uP_Larmor_tramo, a_prop] = TramoAcelPropVariante( v, dt_iter, num_iter, a_prop, a_lab, daprop_diter);
% uP_reac = uP_reac + uP_reac_tramo;
% uP_Larmor = uP_Larmor + uP_Larmor_tramo;
% 
% v_intermedia = v;
% 
% % Desaceleramos, el objetivo es acabar con velocidad 0
% [v, a_lab, uP_reac_tramo, uP_Larmor_tramo, a_prop] = TramoAcelPropVariante( v, dt_iter, num_iter, a_prop, a_lab, -daprop_diter );
% uP_reac = uP_reac + uP_reac_tramo
% uP_Larmor = uP_Larmor + uP_Larmor_tramo
% ~FIN PRUEBA 1

% PRUEBA 2: La aceleración sólo crece --> OK
% dt_iter = 2.5*10^-4;
% num_iter = 2000;
% daprop_diter = [0.8, 0, 0]/num_iter;
% 
% [v, a_lab, uP_reac, uP_Larmor, a_prop] = TramoAcelPropVariante( v, dt_iter, num_iter, a_prop, a_lab, daprop_diter);
% ~FIN PRUEBA 2

% PRUEBA 3: La aceleración sólo crece y luego se mantiene (creciendo la
% velocidad). OK
% dt_iter = 2.5*10^-4;
% num_iter = 2000;
% daprop_diter = [0.8, 0, 0]/num_iter;
% 
% [v, a_lab, uP_reac_tramo, uP_Larmor_tramo, a_prop] = TramoAcelPropVariante( v, dt_iter, num_iter, a_prop, a_lab, daprop_diter);
% uP_reac = uP_reac + uP_reac_tramo;
% uP_Larmor = uP_Larmor + uP_Larmor_tramo;
% 
% v_intermedia = v
% 
% num_iter = 1200;    % Esta fase no será tan larga.
% [v, a_lab, uP_reac_tramo, uP_Larmor_tramo] = TramoAcelPropCte(v, dt_iter, num_iter, a_prop, a_lab);
% uP_reac = uP_reac + uP_reac_tramo;
% uP_Larmor = uP_Larmor + uP_Larmor_tramo;
% ~FIN PRUEBA 3

% PRUEBA 4: aceleración en dos fases en diferentes direcciones... pasable
% :S
dt_iter = 2.5*10^-4;
num_iter = 2000;
daprop_diter = [0.8, 0, 0]/num_iter;

[v, a_lab, uP_reac_tramo, uP_Larmor_tramo, a_prop] = TramoAcelPropVariante( v, dt_iter, num_iter, a_prop, a_lab, daprop_diter);
uP_reac = uP_reac + uP_reac_tramo;
uP_Larmor = uP_Larmor + uP_Larmor_tramo;

v_intermedia = v

daprop_diter = [-1, 0.8, 0]/num_iter;
[v, a_lab, uP_reac_tramo, uP_Larmor_tramo, a_prop] = TramoAcelPropVariante( v, dt_iter, num_iter, a_prop, a_lab, daprop_diter);
uP_reac = uP_reac + uP_reac_tramo;
uP_Larmor = uP_Larmor + uP_Larmor_tramo;
%~PRUEBA 4

v_fin = v
a_prop_fin = a_prop;

PR_fin = PotencialReaccion( v_fin, a_prop_fin )

uP_Tot = uP_reac + uP_Larmor

uP_Tot_menos_PR_fin = uP_Tot - PR_fin

%ratio1 = uP_Tot_menos_PR_fin(1) / uP_reac(1)
%ratio2 = uP_Tot_menos_PR_fin(2) / uP_reac(2)

ratio1 = uP_Tot_menos_PR_fin(1) / uP_Tot(1)
ratio2 = uP_Tot_menos_PR_fin(2) / uP_Tot(2)
