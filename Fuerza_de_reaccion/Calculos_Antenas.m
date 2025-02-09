syms Vo;
syms L;
syms k;
syms x;
syms t;
syms w;

V = Vo*(1-x/L)^2*cos(w*t-k*x);

E = -diff(V,x);

% Densidad de carga
p = diff(E,x);

J = - int(diff(p, t), x)
J2 = - diff(E, t);

Prod = J * E;

%val = int(Prod, x, [0, L])

%w = 1;
%x = 0.5;
%L = 1.5;
%Vo=1;
%t=1.1;
%k=1;

%Jr = subs(J)
%Jr2 = subs(J2)
%Prod_r = subs(Prod)


