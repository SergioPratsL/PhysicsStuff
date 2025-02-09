
c = 299792458;  %m/s

q_e = 1.6 * 10^-19;
m_e = 9.1 * 10^-31;

h = 6.6 * 10 ^-34;

permit = 8.8541 * 10^-12;
permeab = 4 * pi * 10^-7;

q_norm = 1;   % A veces es mejor usar la carga normalizada. 


v_x = 1000;
Ex_ext = 1;

Pot = Ex_ext * v_x * q_norm;

% Recordar que S = E x H = E x b * (1/perm)

Diverg_Poynting = Ex_ext * (q_norm * v_x ) * 2 * 2 / (4*pi);
% El "* 2 * 2" se debe a que el gradiente por un lado es el doble que el
% campo de magnetico y por otro en que hay qye tener en cuenta dos
% dimensiones...


