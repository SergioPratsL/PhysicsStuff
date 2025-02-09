

syms r;
syms a0;

Carga = (4/a0^3) * exp(-2*r / a0) * r^2;

CargaAcum = int(Carga, r)