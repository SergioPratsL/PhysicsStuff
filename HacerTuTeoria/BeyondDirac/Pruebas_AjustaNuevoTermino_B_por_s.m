% La idea de este script es encontrar las constantes que permiten que 
% ε_μijk γ^i ∂^j A^k nos lleva a tener un termino de momento magnético de
% valor -mu * B, siendo que mu es el momento magnético.
% Recordar que el momento magnético del electrón es:
% -ge*(e/2m)*L, siendo ge el factor anómalo que es casi 2 y siendo L=hbar/2

% Debido a las unidades, pense que se podría añadir un término q/mc
% pero habrá que ver, tengo que hacer los cálculos con todas las unidades.

% Lo que aquí se consigue es ver que estas matrices que he propuesto al
% combinarse con un campo magnético da energía en la dirección del spin y
% no la da en las otras direcciones y de bonus, este término hace precesar
% al spin, que en verdad no es poca cosa.

clear;

c_elec = 1.602176 * 10^-19;
m_elec = 9.10938 * 10^-31;
perme = 8.854187 * 10^-12;
c = 299792458;
h = 6.626070 * 10^-34;
h_bar = h / (2*pi);
radio_bohr = 4*pi*perme*h_bar^2 / (m_elec * c_elec^2);
cte_fina = 1/(4*pi*perme) * c_elec^2 /(h_bar*c);

% 08.12 esto no pinta nada aquí pero este fichero se ha vuelto un cajon de
% sastre:
r_espira = h_bar /(m_elec * c * sqrt(3))

[gt, gx, gy, gz] = MatricesGamma();

%g5 = gt * gx * gy * gz;
% Oh no! Un factor "i" diferencia a los dos, gana el MatrizGamma(5)
g5 = MatrizGamma(5);

% Esta es quizá la principal variable con la que juego en este script
spinor_base = [1, 0];               % +Z
%spinor_base = [1, 1] / sqrt(2);     % +X
%spinor_base = [1, i] / sqrt(2);    % +Y

dir_spin = SpinorToVector(spinor_base)

v = [0, 0, 0];
p = fGamma(v) * v;
E = sqrt(1+norm(p)^2);

phi = DiracSpinorPlainWave(p, spinor_base).'

dens_prob = phi' * phi;
phi = phi / sqrt(dens_prob);

phi_gt = phi' * gt;

% Sz = MatrizSpin_4_4(3);
% Sx = MatrizSpin_4_4(1);
% Sy = MatrizSpin_4_4(2);
% contrib = phi_gt * Sz * phi


% Escoger el ordel de las componentes 3 y 4 para dar la combinación qu da B
% positivo, por ejemplo dAx/dy --> +Bz => [x,x,1,2]
if (spinor_base(1) == 1 && spinor_base(2) == 0)
    signoLC = LeviCivita([1,4,2,3]);
elseif (spinor_base(1) == spinor_base(2))
    signoLC = LeviCivita([1,2,3,4]);
else
    signoLC = LeviCivita([1,3,4,2]);
end

% La idea de este cálculo es estudiar la contribución de un 
% campo magnético sober el Hamiltoniano debido al nuevo término.
% contrib2t en verdad no pinta nada. 
contrib2z = - phi_gt * g5 * gz * phi;
contrib2x = - phi_gt * g5 * gx * phi;
contrib2y = - phi_gt * g5 * gy * phi;
contrib2t = - phi_gt * g5 * gt * phi;


% Si el spin es en X, la única contribución que no es cero es la
% contribución en X, si es en Y es la Y, si es en Z, la Z. La t siempre
% es cero y además todos los valores son +1i
% ¡RESULTADO IMPORTANTE!
contribuciones = signoLC * [contrib2t, contrib2x, contrib2y, contrib2z].';

% matriz_z = 1i * (g5 * gz);
% matriz_x = 1i * g5 * gx;
% matriz_y = g5 * gy;
% matriz_t = g5 * gt;
% 
%  xxx = - 1i * gt * matriz_y;
% % St_obtenido_por_induccion = - 1i * gt * g5 * gt:
% 
% % Anticonmutan
% Gamma_5_Gamma_t_conmut = gt * g5 + g5 * gt;

if (spinor_base(1) == 1 && spinor_base(2) == 0)
    g_dir = gz;
elseif (spinor_base(1) == spinor_base(2))
    g_dir = gx;
else
    g_dir = gy;
end


% La hora de la verdad, esto cambia los autovectores... supongamos que el
% campo es dxAy --> +Bz:
%termino_H = -1i * gt * g5 * gz * phi
%termino_H = gt * g5 * gz * phi;
termino_H = gt * g5 * g_dir * phi;

% El termino_H es proporcial a Phi si phi es en Z

dir_cambio_spin = SpinorToVector(termino_H(1:2).')

% 12.10. Cambie g5 y no revisé esto...
% Una vez le quitas el "i" las contribuciones hacen que precese sobre el
% plano XY, la dir_spin sigue siendo  en sentido opuesto (para X e Y)
% pero el loop siguiente muestra que el bicho precesa :).

% phi_actual = phi;
% factor_dif = 1/1000;
% %matriz_nuevo_termino = -1i * gt * g5 * gz;
% matriz_nuevo_termino = gt * g5 * gz;
% 
% n = 1;
% while n <= 1000
%     termino_H = matriz_nuevo_termino * phi_actual;
%     
%     phi_actual = phi_actual + termino_H / factor_dif;
%     phi_actual = phi_actual / norm(phi_actual);
%     
%     n = n + 1;
% end
% 
% phi_final = phi_actual
% dir_spin_final = SpinorToVector(phi_final(1:2).')


% Contribución al momento causada por el campo magnético
% en las tres direcciones p_dir_prt = -i phi
p_dir_prt = g_dir * g5 * gt * phi


