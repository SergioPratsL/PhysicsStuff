function [factor_atenuacion, factor_inicial, ds_iteracion] = CalculaFactorAtenuacion(dt_iteracion, Sigma)

ds_iteracion = dt_iteracion / Sigma;

% Tener en cuuenta que el periodo en el SRI propio es más corto, tened en
% cuenta que dt_iteracion es negativo
factor_atenuacion = exp( (3/2)*ds_iteracion ); 

factor_inicial = 1 - factor_atenuacion;

end

