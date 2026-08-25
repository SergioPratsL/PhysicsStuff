function [v_fin, acel_fin, uP_reac, uP_Larmor] = TramoAcelPropCte( v_ini, dt_iter, num_iter, a_prop, a_lab_ini )

% Parametros de entrada:
  % v_ini: la velocidad de inicial de la partícula.
  % dt_iter: el tiempo de laboratorio que queremos que pase entre dos
  % iteraciones
  % num_iter: el número de iteraciones que pasa aquí
  % a_prop: el ritmo de aceleración propia que es constante.
  % a_lab_ini es la aceleración del laboratorio inicial --> necesario

% Parámetros de salida:
  % v_fin: la velocidad de la partícula tras haber estado acelerando.
  % acel_fin: la aceleración de laboratorio de la partícula al acabar.
  % desde el laboratorio..
  % uP_react: es el 4 vector de energía y momento causado por la fuerza de
  % reacción
  % uP_Larmor: es el 4 vector de energía y momento entregada al campo
  % radiado
  
  uP_reac = [0,0,0,0];
  uP_Larmor = [0,0,0,0];
  
  v = v_ini;
  a = a_lab_ini;
  
  dp = a_prop * dt_iter;
  
  n = 1;
  while n <= num_iter
     
     % De esta forma no hacen falta las variables "old" porque la misma 
     %variable se pasa y al recuperarse ya se ha vuelto vieja (servirá 
     % para la siguiente iteración siendo entonces vieja) 
     [v, a, uP_reac, uP_Larmor] = AccionesPotencialReaccion(v, dt_iter, a_prop, a, uP_reac, uP_Larmor);
     
     n = n + 1;
  end
  
  v_fin = v;
  acel_fin = a;
