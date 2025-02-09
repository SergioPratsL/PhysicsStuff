function [E_uv, B_uv] = fCrea_Matriz_Derivadas( Cofs )
% Esta función obtiene la matriz con las derivadas de E y B en cada una de
% las coordenadas, la "u" representa la dirección de las derivadas y la "v"
% la componente del campo por lo que si v=0, siempre tenemos 0 (no existe
% Et o Bt). Los coeficientes siguen un criterio arbitrario bajo el cual he
% elegido 16 componentes libres y luego otras son calculadas a partir de
% estas (básicamente siempre son calculadas las del campo magnético).

% Componentes del vector:
% [dxEx, dxEy, dtEx, dtEy, dtEz, dxEy, dxEz, dyEx
% dyEz, dzEx, dzEy, dxBx, dyBy, dxBy, dyBz, dzBx]

% Pendiente ver cómo se trata la carga y la corriente (espero que no sean necesarias :P).

E = zeros(4);
B = zeros(4);

E(2,2) = Cofs(1);
E(3,3) = Cofs(2);
E(4,4) = -(E(2,2) + E(3,3));        % Dens carga = 0

E(1,2) = Cofs(3);
E(1,3) = Cofs(4);
E(1,4) = Cofs(5);

E(2,3) = Cofs(6);
E(2,4) = Cofs(7);
E(3,2) = Cofs(8);

E(3,4) = Cofs(9);
E(4,2) = Cofs(10);
E(4,3) = Cofs(11);

B(2,2) = Cofs(12);
B(3,3) = Cofs(13);
B(4,4) = -(B(2,2) + B(3,3));

B(1,2) = - (E(3,4)-E(4,3));
B(1,3) = - (E(4,2)-E(2,4));
B(1,4) = - (E(2,3)-E(3,2));

B(2,3) = Cofs(14);
B(3,4) = Cofs(15);
B(4,2) = Cofs(16);

B(3,2) = B(2,3) - E(1,4);   % J=0
B(2,4) = B(4,2) - E(1,3);   % J=0
B(4,3) = B(3,4) - E(1,2);   % J=0


E_uv = E;
B_uv = B;



end

