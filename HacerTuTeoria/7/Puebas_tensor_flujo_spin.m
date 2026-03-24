
c_elec = 1.602176 * 10^-19;
m_elec = 9.10938 * 10^-31;
perme = 8.854187 * 10^-12;
c = 299792458;
h = 6.626070 * 10^-34;
h_bar = h / (2*pi);
cte_fina = 1/(4*pi*perme) * c_elec^2 /(h_bar*c);
permeabilidad = 2 * cte_fina / c_elec^2 * (h/c);
radio_bohr = 4*pi*perme*h_bar^2 / (m_elec*c_elec^2);
Z = 1;

factor_gamma_wiki = sqrt(1-cte_fina^2);    


[gt, gx, gy, gz] = MatricesGamma();

%x = MatrizSpin_4_4(1)
%Sx_2 = i/2 * (gy*gz - gz*gy)

%Sy = MatrizSpin_4_4(2)
%Sy_2 = i/2 * (gz*gx - gx*gz)

%Sz = MatrizSpin_4_4(3)
%Sz_2 = i/2 * (gx*gy - gy*gx)

Phi = [1;0;0;0];

% Quiero mirar primero el T10-T01 en dirección X, luego T20-T02 en Y, y
% luego T30-T03 en Z.

% En cada dirección el 'a' vale -.5 y el 'b' + .5.
val_xa = TensorFlujoSpin(Phi, 1, 0, 1);
val_xb = TensorFlujoSpin(Phi, 0, 1, 1);

val_ya = TensorFlujoSpin(Phi, 2, 0, 2);
val_yb = TensorFlujoSpin(Phi, 0, 2, 2);

val_za = TensorFlujoSpin(Phi, 3, 0, 3);
val_zb = TensorFlujoSpin(Phi, 0, 3, 3);


% Valor que aplica para el problema de hidrógeno, sea la direcciónde r = x
% Necesito la verdadera forma de la onda

spinor_base = [1, 0];

dir = [1,0,0];

% Con la base ya basta para este caso
bispinor_base = MontaBispinorBase(factor_gamma_wiki, cte_fina, dir, Z).';
bispinor_base = bispinor_base / norm(bispinor_base);

% Real, 0 e imaginario... Podría funcionar...
val_y = TensorFlujoSpin(bispinor_base, 0, 2, 1);
val_z = TensorFlujoSpin(bispinor_base, 0, 3, 1);
val_x = TensorFlujoSpin(bispinor_base, 0, 1, 1);


dir = [0,1,0];
bispinor_base = MontaBispinorBase(factor_gamma_wiki, cte_fina, dir, Z).';
bispinor_base = bispinor_base / norm(bispinor_base);

% Imaginario, 0 y real... ok...
val_y = TensorFlujoSpin(bispinor_base, 0, 2, 2)
val_z = TensorFlujoSpin(bispinor_base, 0, 3, 2)
val_x = TensorFlujoSpin(bispinor_base, 0, 1, 2)


function bispinor_base = MontaBispinorBase(factor_gamma_wiki, cte_fina, dir, Z)
    spinor_grande = (1+factor_gamma_wiki) * [1,0];
    
    spinor_pequeno = -1i * Z * cte_fina * [dir(3), (dir(1)+1i*dir(2))];

    bispinor_base = [spinor_grande, spinor_pequeno];
end


