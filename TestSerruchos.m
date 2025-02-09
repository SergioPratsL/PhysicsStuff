% Matlab para ver si un gas puede hacer fuerza contra los famosos
% serrunchos pratianos, que desgraciadamente solo existen en mi mente

% Las distancias dadas en verdad son en nanometros
clear

D = 100;
H = 20;
%alfa = pi / 6;
%alfa = pi / 4;
%alfa = pi / 3;
alfa = pi / 12;

% Longitud del serrucho
tan_alfa = tan(alfa)
L = H / tan_alfa;

cos_alfa = cos(alfa);
sen_alfa = sin(alfa);

numRho = 100;
%numRho = 200;
%numPsi = 50;
%numPsi = 100;
%numPsi = 400;
%numPsi = 1000;
numPsi = 100;
%num_l = 10;        
%num_l = 50;        
num_l = 20;

l_unit =  L / (num_l + 1);
% Rho se integra entre 0 y pi /2
Rho_unit = (pi/2) / (numRho + 1);
% Psi se integra entre 0 y pi. (No nos importa si la z es positiva o
% negativa
Psi_unit = pi / (numPsi + 1);

l_ini = l_unit / 2;
Rho_ini = - Rho_unit / 2;       % Va hacia abajo
Psi_ini = Psi_unit / 2;

Presion_x_acum = 0;
contador_bolas = 0;
contador_solid = 0;

iter_l = 0;
while iter_l < num_l
    l = l_ini + l_unit * iter_l;
    pos = [l, D];
    
    iter_Rho = 0;
    while iter_Rho < numRho
        Rho = Rho_ini - Rho_unit * iter_Rho;    % El menos por ir hacia abajo!
        
% Factor de normalizacion, me como pi's y otras mierdas, lo que importa es
% este seno
        ang_solid = cos(Rho);
    
        iter_Psi = 0;
        while iter_Psi < numPsi
            Psi = Psi_ini + Psi_unit * iter_Psi;
            
% Salvaguarda.            
            if Psi == pi / 2
                Psi = pi / 2 - 0.0001;
            end
            
            
%            p_x = ReboteSerruchos( pos, Rho, Psi, 0, H, L, alfa, tan_alfa );
            if iter_l == 0 && iter_Rho == 12 && iter_Psi == 26
                dsadaS = 59;
            end

            p_x = ReboteSerruchos( pos, Rho, Psi, H, L, tan_alfa, cos_alfa , sen_alfa );            
            
            Presion_x_acum = Presion_x_acum + p_x * ang_solid;
            
            contador_bolas = contador_bolas + 1;
            contador_solid = contador_solid + ang_solid;
            iter_Psi = iter_Psi + 1;
        end
        iter_Rho = iter_Rho + 1;
    end
    iter_l = iter_l + 1;
end


contador_solid = contador_solid
Presion_x_acum = Presion_x_acum 

ratio = Presion_x_acum  / contador_solid


