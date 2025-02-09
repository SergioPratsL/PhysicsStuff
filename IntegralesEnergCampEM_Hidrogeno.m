% Comento conforme no necesito más el uso de los campos
clear

syms r a L;

Onda_quad_rad_1 = 4*pi*r^2*exp(-2*r/a)/ (a^3 * pi);

%int_quad = int( Onda_quad_rad_1, r, 0, 1000 )

dens_e_pot_rad_1 = Onda_quad_rad_1 / r;

E_pot_1 = int(dens_e_pot_rad_1, r, 0, 1000*a)

 %Carga_Int_1 = int(Onda_quad_rad_1, r, 0, L)
% Resultado:   (Esta la clavé!!)
% Carga_Int_1 = 1 - (exp(-(2*L)/a)*(2*L^2 + 2*L*a + a^2))/a^2
    
%Campo_E_1 = Carga_Int_1 / L^2

%Prod_Campos_1 = - 2 * Campo_E_1 * (1 / L^2);  
%Prod_Campos_rad_1 = Prod_Campos_1 * 4 * pi * L^2;

%Ener_tot_1 = int( Prod_Campos_rad_1, L, 0, 1000 )

% Ener_elect_E_1 = 8*pi*(limit(- exp(-(2*L)/a)/a - (exp(-(2*L)/a) - 1)/L, L, Inf) - 1/a)

%Dens_Energ_E_rad_1 = Campo_E_1^2 * 4 * pi * L^2;
%Ener_elect_E_1 = int( Dens_Energ_E_rad_1, L, 0, 20000*a )












