function [v_fin, acel_fin, uP_reac, uP_Larmor, a_prop_fin] = TramoAcelPropVariante( v_ini, dt_iter, num_iter, a_prop, a_lab_ini, daprop_diter )

% Parametros de entrada:
  % v_ini: la velocidad de inicial de la partícula.
  % dt_iter: el tiempo de laboratorio que queremos que pase entre dos
  % iteraciones
  % num_iter: el número de iteraciones que pasa aquí
  % a_prop: el ritmo de aceleración propia que es constante.
  % a_lab_ini es la aceleración del laboratorio inicial --> necesario
  % daprop_diter: Cuanto cambia la aceleracion propia por iteracion,
  % úsalo con sabiduría

% Parámetros de salida:
  % v_fin: la velocidad de la partícula tras haber estado acelerando.
  % acel_fin: la aceleración de laboratorio de la partícula al acabar.
  % desde el laboratorio..
  % uP_react: es el 4 vector de energía y momento causado por la fuerza de
  % reacción
  % uP_Larmor: es el 4 vector de energía y momento entregada al campo
  % radiado
  % a_prop_fin: la aceleración propia que queda al final.
  
  if(norm(daprop_diter) == 0)
      disp('Si daprop_diter = 0, usa TrampoAcelPropCte en vez de esta funcion!')
  end
  
  uP_reac = [0,0,0,0];
  uP_Larmor = [0,0,0,0];
  
  v = v_ini;
  a = a_lab_ini;
  
  n = 1;
  while n <= num_iter
      
     if n == 1500
         n = 1500;
     end

     a_prop = a_prop + daprop_diter;
     [v, a, uP_reac, uP_Larmor] = AccionesPotencialReaccion(v, dt_iter, a_prop, a, uP_reac, uP_Larmor);

     n = n + 1;
  end
  
  v_fin = v;
  acel_fin = a;
  a_prop_fin = a_prop;