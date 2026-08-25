function [v, a, uP_reac, uP_Larmor] = AccionesPotencialReactivo( v_old, dt, a_prop, a_old, uP_reac_old, uP_Larmor_old )

% Esta función sirve para abstraer el bucle de Acel_Prop_Cte y de otras
% funciones similares que llegarán

dp = dt * a_prop;

[v, Sigma] = ModificaVelocidad(v_old, dp);
a = (v - v_old) / dt;
da_dt = (a - a_old) / dt;
     
Fr = FuerzaReaccionCovariante(dt, v, a, da_dt, Sigma);
duP_reac = [dot(v, Fr), Fr] * dt;
uP_reac = uP_reac_old + duP_reac;
     
%duP_Larmor = LarmorCovarianteMomentum(v, a_prop, dt, Sigma);
duP_Larmor = LarmorCovarianteMomentum2(v, a_prop, dt);
uP_Larmor = uP_Larmor_old + duP_Larmor;
