
function Rr = DRetardada( R, v )

% Funcion para obtener la distancia retardada a partir de la distancia
% actual

% R es la distancia actual desde el emisor al punto a medir
% v es la velocidad con la que se mueve el emisor expresada en partes por
% unidad de la velocidad de la luz

% Esta version asumira que la velocida no tiene que ser en el eje X...
% puedo dar mas que eso!

% Hay 3 ecuaciones lineales mas una cuadratica, esto requiere un bolzano
% o bien definir una base vectorial


if norm(v) == 0
    Rr = R
    return
end

speed = norm(v);
v_norm = v / speed;


% Ahora elijo la segunda direccion mediante un criterio arbitrario pero
% fiable

x = [1, 0, 0];

v_aux1 = cross( v_norm, x );

% Caso en que la velocidad es  en X
if norm(v_aux1) == 0
    v_aux1 = [0, 1, 0];
else
    v_aux1 = v_aux1 / norm(v_aux1);
end


v_aux2 = cross( v_norm, v_aux1 );

% Ahora a descomponer el vector en su nueva base
val_v_norm = dot( R, v_norm );
val_v_aux1 = dot( R, v_aux1 );
val_v_aux2 = dot( R, v_aux2 );

Sigma = 1/sqrt(1-speed^2);

raiz_chunga = sqrt( Sigma^2 * val_v_norm^2 + val_v_aux1^2 + val_v_aux2^2 );

% Ojo al uso del signo mas en la raiz!!
val_v_norm_mod = Sigma^2 * val_v_norm + speed * Sigma * raiz_chunga;


% Ahora a componer el vector de salida!!

Rr = val_v_norm_mod * v_norm + val_v_aux1 * v_aux1 + val_v_aux2 * v_aux2;

% Parte de verificacion:

dist_Rr = norm(Rr);

Rr_verif = R + dist_Rr * v;

dif_verif = Rr - Rr_verif;






