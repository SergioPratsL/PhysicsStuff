clear;

% El objetivo de este scrip es confirmar mis sospechas que las 6 igualdades
% de segundo grado que sirven para equilibrar los flujos del tensor E-E de
% Dirac en verdad no son más que identidades para una elección de Phi y
% dPhi/dx dada por el Hamiltoniano. De confirmarse, debo abandonar toda
% esperanza puesta en ellas.

% Estas 6 ecuaciones son las que llamaba "jinetes del apocalipsis de
% Dirac", nótens que esto aplica a un escenario como el QHO donde H
% conmutas con 2 direcciones y por tanto sólo queda una dirección a
% determinar (X). En este escenario dada una energía y dada X, phi
% determina a dphi/dx o viceversa, entonces podemos elegir valores de X, E
% varios y valores de phi, ello nos fuerza a elegir dphi/dx, entonces con
% todos los valores vemos si las ecuaciones se cumplen o no. Si siempre se
% cumplen, son identidades que emanan del Hamiltoniano y por tanto no me
% pueden ayudar en nada, si no se cumplen sí que podrían ayudar y habría
% que resolverlas.

% Cabe destacar que en caso que las identidades se cumplan, debo echar un
% nuevo vistazo al azote del Diario 6 (y finalmente de todos), al
% QHO_Dirac_test.m para ver qué hice y que pasa, porque significaría que
% con H y con E (E a determinar en cada prueba), tienes todo lo que hace
% falta.

% RIESGO: que me equivoque y por culpa de ello las identidades no se
% cumplan


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


% Hay que ir probando muchos valores y ver qué sale.

% Valor en el que quiero evaluar la onda
% x = 0.25 * radio_bohr;
% %phi = [1; 0; 0.1*i; 0];
% %phi = [1; 0; 0.1; 0];
% %phi = [1; 0.3*i; 0.2; -0.15*i];
% phi = [0.9; i; 0.4i; 0.2*i];
% E = mc2 + nivel_base;

% x = 0;
% %phi = [1; 0; 0.1*i; 0];
% %phi = [1; 0; 0.1; 0];
% %phi = [1; 0.3*i; 0.2; -0.15*i];
% phi = [0.9; i; 0.4i; 0.2*i];
% E = mc2 + 2.35 * nivel_base;

x = 3.15 * radio_bohr;
%phi = [1; 0; 0.1*i; 0];
%phi = [1; 0; 0.1; 0];
%phi = [1; 0.3*i; 0.2; -0.15*i];
phi = [0.9; i; 0.4i; 0.2*i];
E = mc2 + 8 * nivel_base;

phi = phi / norm(phi);
E_pot = -1/2 * k * x^2;

E_can = E - E_pot;

% Para determinar las derivadas del spinor "B" en función de las
% componentes del spinor
cA = -1i/hc*(E - mc2 - E_pot);
cB = -1i/hc*(E + mc2 - E_pot);

% Determinar la derivada para el QHO en función de estos parámetros.
%%%dPhi_dx = [cB * phi(3:4); cA * phi(1:2)];
dPhi_dx = [cB*phi(4); cB*phi(3); cA*phi(2); cA*phi(1)];

% ¡¡¡ En este contexto dPhi_dy = dPhi_dz = 0 !!!


% Se usa en ecuaciones.

%M_xtx = (1i/2) * gx * (gx*gt - gt*gx);
%M_ytx = (1i/2) * gx * (gy*gt - gt*gy);
%M_ztx = (1i/2) * gx * (gz*gt - gt*gz);

% Como dijo Pique, el gt "Se queda".
M_xtx = (1i/2) * gt * gx * (gx*gt - gt*gx);
M_ytx = (1i/2) * gt * gx * (gy*gt - gt*gy);
M_ztx = (1i/2) * gt * gx * (gz*gt - gt*gz);

M_xyx = (1i/2) * gt * gx * (gx*gy - gy*gx);
M_xzx = (1i/2) * gt * gx * (gx*gz - gz*gx);
M_yzx = (1i/2) * gt * gx * (gy*gz - gz*gy);

ParteIzda_1 = E_can/c * phi' * ax * phi - i * h_bar * phi' * dPhi_dx;
ParteDerecha_1 = 1/2 * h_bar * (dPhi_dx' * M_xtx * phi + phi' * M_xtx * dPhi_dx);
Resultado_1 = ParteIzda_1 - ParteDerecha_1


ParteIzda_2 = E_can/c * phi' * ay * phi;
ParteDerecha_2 = 1/2 * h_bar * (dPhi_dx' * M_ytx * phi + phi' * M_ytx * dPhi_dx);
Resultado_2 = ParteIzda_2 - ParteDerecha_2


ParteIzda_3 = E_can/c * phi' * az * phi;
ParteDerecha_3 = 1/2 * h_bar * (dPhi_dx' * M_ztx * phi + phi' * M_ztx * dPhi_dx);
Resultado_3 = ParteIzda_3 - ParteDerecha_3


% No quito unidades para que todo esté en el mismo órden de magnitud.
% Que estos tengan una "i" para el lado derecho y los otros dos es extraño,
% pero debe ser causadoo por ser dos componentes espaciales en vez de una.
ParteIzda_4 = h_bar * phi' * ay * dPhi_dx;
ParteDerecha_4 = i * h_bar / 2 * (dPhi_dx' * M_xyx * phi + phi' * M_xyx * dPhi_dx);
Resultado_4 = ParteIzda_4 - ParteDerecha_4

ParteIzda_5 = h_bar * phi' * az * dPhi_dx;
ParteDerecha_5 = i * h_bar / 2 * (dPhi_dx' * M_xzx * phi + phi' * M_xzx * dPhi_dx);
Resultado_5 = ParteIzda_5 - ParteDerecha_5

ParteIzda_6 = 0;
ParteDerecha_6 = i * h_bar / 2 * (dPhi_dx' * M_yzx * phi + phi' * M_yzx * dPhi_dx);
Resultado_6 = ParteIzda_6 - ParteDerecha_6



% Parece que 3 de los seis jinetes bastaron para apalizar a Sergio Plato,
% pues los jinetes quizá son necesarios.... Parece que no.


