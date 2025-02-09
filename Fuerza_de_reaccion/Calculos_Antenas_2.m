syms w;
syms t;
syms k;
syms x;
syms L;
syms E0;

% Este vector representa: [cos(wt), cos(wt-kL), sen(wt), sen(wt-kL)]
Values = [-1,0,0,0];

% Los generadores son [sen(wt-kx), (x/L)*cos(wt-kx), (x/L)*sen(wt-kx), cos(wt-kx)]

Cofs1 = [-1/k,-1/(L*k^2),0,0];
Cofs2 = [1/k,1/(L*k^2),1/k,0];
Cofs3 = [0,0,-1/(L*k^2),1/k];
Cofs4 = [0,1/k,1/(L*k^2),-1/k];

Cofs = [Cofs1;Cofs2;Cofs3;Cofs4];


%Result = linsolve(Cofs, Values')

E = E0 * (1 - x / L) * sin(w*t-k*x);

Int_E = int(E, x)

V = simplify(subs(Int_E, x, L) - subs(Int_E, x, 0) );

J = -diff(E, t);

J0 = subs(J, x, 0);

Prod_0 = V2 * J0;

Int_Prod_0 = int(Prod_0, t);

Int_Prod_0_Period = simplify(subs(Int_Prod_0, t, (2*pi()/w)) - subs(Int_Prod_0, t, 0) );

p = diff(E, x)


