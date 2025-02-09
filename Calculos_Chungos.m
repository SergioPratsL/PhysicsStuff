syms r;

syms a;

f = r^2 * exp(-r/a);

Fr = int(f,r)

% Fr = -a^3*exp(-r/a)*((2*r)/a + r^2/a^2 + 2)