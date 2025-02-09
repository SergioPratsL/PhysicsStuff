syms x
syms y
syms z

f = sqrt(z^2/(x^2+y^2+z^2))   % (sqrt(x^2/(x^2+y^2))+1*i*sqrt(y^2/(x^2+y^2)));

df_dx = diff(f, x);
df_dx2 = diff(df_dx, x);

df_dy = diff(f, y);
df_dy2 = diff(df_dy, y);

df_dz = diff(f, z);
df_dz2 = diff(df_dz, z);

df_dxdz = diff(df_dx, z);

x = 1/sqrt(2)*sqrt(3)/2;
y = 1/sqrt(2)*1/2;
z = 1/sqrt(2);

val_dx2 = subs(df_dx2)
val_dy2 = subs(df_dy2)
val_dz2 = subs(df_dz2)
val_dxdz = subs(df_dxdz)

val_tot = val_dx2 + 2 * val_dy2 + val_dz2 + 2 * val_dxdz;

% Factor que hay que aplicar al final.
val_tot = - val_tot / 2