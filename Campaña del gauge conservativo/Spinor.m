function u = Spinor( p, Signo )
% Se asume que la masa esta normalizada
% para seguir nomenclatura de Feynman
% Basado en formulas paginas 58-60 del libro
% Es una simplificacion que no contempla pz <> 0

mom = dot(p,p);

E = sqrt(1 + mom);

F = E + 1;

p_mas = p(1) + i* p(2);
 
p_menos = p(1) - i* p(2);


if Signo == '+'
    u = [F, 0, 0, p_mas];
elseif Signo == '-'
    u = [0, F, p_menos, 0];
end

% Dejarlo en vertical
u = u';
u = conj(u);



%p = [0.3, 0.1, 0]
% u = Spinor ( p, '+')
%  v =  u';
% jx = v * M_x * u

% Haciendo esto que diria que esta bien me sale que la onda plana tiene
% corriente
 
% Si quitase M_t (me la olvide) y conjugase v, entonces da cero...

% v = conj(v);
% jx = v * M_t * M_x * u   (de nuevo)

% jy = v * ( M_t * ( M_y * u ) )


% Mas misticismo.

% La ratio entre jx y jy es 3, pero si cambio el spin es -3...


