clear

syms r a L;

Onda_quad_rad_2 = 4*pi*r^2*(2-r/a)^2 * exp(-r/a)/ a^3 / (32*pi);

%int_quad = int( Onda_quad_rad_2, r, 0, 1000 )

dens_e_pot_rad_2 = Onda_quad_rad_2 / r;

E_pot_2 = int(dens_e_pot_rad_2, r, 0, 1000*a)


%Carga_Int_2 = int(Onda_quad_rad_2, r, 0, L)
% Resultado:   
% Carga_Int_2 = 1 - (exp(-L/a)*(L^4 + 4*L^2*a^2 + 8*L*a^3 + 8*a^4))/(8*a^4)

%Campo_E_2 = Carga_Int_2 / L^2

%Prod_Campos_2 = - 2 * Campo_E_2 * (1 / L^2);  
%Prod_Campos_rad_2 = Prod_Campos_2 * 4 * pi * L^2;

%Ener_tot_2 = int( Prod_Campos_rad_2, L, 0, 1000 )
% - (2*pi)/a - (pi*exp(-1000/a)*(2000000*a - 8*a^3*exp(1000/a) + 6000*a^2 + 8*a^3 + 1000000000))/(1000*a^3)


%Dens_Energ_elec_E_rad_2 = Campo_E_2^2 * 4 * pi * L^2;

%Ener_E_2 = int( Dens_Energ_E_rad_2, L, 0, 50000000*a )



