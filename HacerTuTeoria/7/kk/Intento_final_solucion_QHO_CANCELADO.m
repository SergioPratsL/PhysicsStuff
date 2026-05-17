%% CANCELADO TRAS DEMOSTRARSE QUE LOS SEIS JINETES DEL APOCALIPSIS
%% DE DIRAC NO ERAN MÁS QUE 6 IMPOSTORES!

clear;

% El objetivo es poder demostrar que el QHO no tiene solución en lugar de
% simplemente decir "no sé".

q = 1.602176 * 10^-19;
m_e = 9.10938 * 10^-31;
perme = 8.854187 * 10^-12;
c = 299792458;
h_bar = 6.626070 * 10^-34 / (2*pi);
cte_fina = 1/(4*pi*perme) * q^2 /(h_bar*c);
%permeabilidad = 2 * cte_fina / q^2 * (h/c);
radio_bohr = 4*pi*perme*h_bar^2 / (m_e*q^2);

mc2 = m_e * c^2;
hc = h_bar * c;

[gt, gx, gy, gz] = MatricesGamma();
g5 = MatrizGamma(5);

Sigma_xt = gx*gt - gt*gx;


[at, ax, ay, az] = MatricesAlfa();

%[sx, sy, sz] = MatricesPauli();

Sx = MatrizSpin_4_4(1);

% k se usa para determinar el potencial
% Los valores deben ser negativos porque en un potencial atractivo son
% negativos, cuanto más lejos, más restan.
k = -1000;

% para que el factor potencial sea 1/2*w^2*x^2
w_wiki = sqrt(k/m_e);
nivel_base = h_bar * w_wiki / 2;


% dx tiene que ser suficientemente pequeña para que en su transcurso la
% derivada no cambie demasiado y así vemos si los valores obtenidos en x y
% x+dx son coherentes o no.
dx = radio_bohr / 1000;

% Energia a la que evaluo la onad (el parámetro que anda suelto).
E = mc2 + nivel_base;



% Valor en el que quiero evaluar la onda
x = 0;

E_pot = -1/2 * k * x^2;

% No sirve para nada, sólo es orientativo
V = k * x^2 /q;


% Para determinar las derivadas del spinor "B" en función de las
% componentes del spinor
cA = -1i/hc*(E - mc2 - E_pot);
cB = -1i/hc*(E + mc2 - E_pot);


% By Gemini:

% % Ver qué puede hacer el Matlab.
% phi = sym('phi', [4 1]);
% phi_dag = conj(phi.');
% 
% Jx = phi_dag * ax * phi;
% 
% % The x-current density expression
% Jx = phi_dag * ax * phi;
% 
% % Use 'expand' to see every individual quadratic term
% Jx_expanded = expand(Jx);
% 
% disp('Expanded x-current terms:');
% pretty(Jx_expanded)



% Ver cómo hago los vectores que se montan a partir de phi

syms w x y z
phi = sym([w;x;y;z]);
phi_camb = sym([cA*w;cA*x;cB*y;cB*z]);

phi_conj = conj(phi.');
phi_camb_conj = conj(phi_camb.');

% Matrices para la divergencia del flujo de Suop...
M_xtx = (1i/2) * gx * (gx*gt - gt*gx);

% Esta sólo es la primera parte de la ecuación...
Ec1 = mc2 * phi_conj * ax * phi;




