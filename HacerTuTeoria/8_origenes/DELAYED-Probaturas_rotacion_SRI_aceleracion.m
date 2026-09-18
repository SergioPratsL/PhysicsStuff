% Este Matlab se basa en asumir que la rotación de los sistemas 
% de referencia acelerados en base a la precesión W-T con el SRIC es algo
% verdadero (pretendo investigar que es algo que funciona bien)

% El objetivo de estas pruebas es sacar la aceleración propia a partir de
% la aceleración vista por un fulano B que no es el SRI pero que conoce al
% SRIC y por tanto es capaz de hacer ciertas cábalas. Sea "A" la particula
% que acelera, B el observador y C el SRIC, empezará sambiendo vBA, vBC y
% dv' que es la aceleración en dt' vista por B (nótese que la comilla se
% asigna a B, dado que es un extraño.


% OJO: ¡A lo mejor permitir cualquier aceleración para un vCA dado no es
% válido! ¡En tal caso me estaría saliendo del dominio de validez de estas
% hipótesis y me metería una hostia!