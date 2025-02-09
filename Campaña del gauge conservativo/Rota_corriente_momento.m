function Js_o = Rota_corriente_momento( Js_i, ang)

% Solo admite rotaciones en el plano XY
% ang en radianes!!

% Valor intermedio
Js_m = Js_i;

% Primero rotar las coordenadas de cada vector

% Nota: no se si esta formula tiene el signo menos mal pero en tal caso
% seria como rotar con la direccion contraria

Js_m(1,1) = cos(ang) * Js_i(1,1) - sin(ang) * Js_i(2,1);
Js_m(2,1) = cos(ang) * Js_i(2,1) + sin(ang) * Js_i(1,1);

Js_m(1,2) = cos(ang) * Js_i(1,2) - sin(ang) * Js_i(2,2);
Js_m(2,2) = cos(ang) * Js_i(2,2) + sin(ang) * Js_i(1,2);

Js_m(1,3) = cos(ang) * Js_i(1,3) - sin(ang) * Js_i(2,3);
Js_m(2,3) = cos(ang) * Js_i(2,3) + sin(ang) * Js_i(1,3);


% Ahora rotar lo que es la corriente de X y la e Y:

Js_o = Js_m;

Js_o(:,1) = cos(ang) * Js_m(:,1) - sin(ang) * Js_m(:,2);
Js_o(:,2) = cos(ang) * Js_m(:,2) + sin(ang) * Js_m(:,1);



