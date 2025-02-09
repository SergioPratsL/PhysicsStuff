% Este script pretende evaluar el campo radiado por una superficie esférica
% de radio 1 sobre un punto que está en el eje X a una distancia arbitraria
% y ante una aceleración en cualquier dirección.

% Se entiende que la superficie esférica está en reposo y acelera toda ella
% de manera uniforme.

% El objetivo de este script es ver si a distancias pequeñas el campo
% radiado se comporta como si la carga acelerada estuviese en el centro de
% la esfera o no.

Aceleracion = [0, 1, 0];
%Aceleracion = [0.5, sqrt(3)/2, 0];
%Aceleracion = [-.005, 1, 0] / norm([.01, 1, 0]);
%Aceleracion = [0, 1, -.005] / norm([.01, 1, 0]);
Velocidad = [0, 0, 0],
RadioEsfera = 1;

%DistanciaX = 0.5
%DistanciaX = 0.1;
%DistanciaX = -0.005;
DistanciaX = 0.0001;
%DistanciaX = 999;      % Converge con carga puntual como era de esperar.
%DistanciaX = -0.2;      % Dentro de la esfera en  x=0.8
%DistanciaX = -0.5;
%DistanciaX = -0.8;
PosicionEvaluada = [1+DistanciaX, 0, 0];

% Cuando procesemos puntos por debajo de esta distancia (si el objeto está
% muy cerca) reduciremos la resolución para enfatizar las contribuciones d
% de los puntos más cercanos
DistanciaSiguienteReduccionIni = 0.25;
DistanciaSiguienteReduccion = DistanciaSiguienteReduccionIni;

% En verdad el Thetha hace cortes a X constante, no a Y constante!!
NumIterThethaIni = 200;
%NumIterThethaIni = 400;
NumIterThetha = NumIterThethaIni;
% Rho reparte la fase entre Y y Z para cada plano, la Z juega el papel de
% la X!!
NumIterRho = 180;
%NumIterRho = 360;
NumIterRhoIni = NumIterRho;

CofReductorDistancia = 4;
CofIncremNumIterThetha = 2;
CofIncremNumIterRho = 2;

ResolRho = (2*pi) / NumIterRho;
ResolThetha = pi / NumIterThetha;

num_iteraciones_hechas = 0;
superficie_acumulada = 0;
E_rad_tot = [0,0,0];
B_rad_tot = [0,0,0];

% Empezamos por Theta cercano a 180º
ThethaInicial = pi - ResolThetha/2;

Thetha = ThethaInicial;
while Thetha > 0
    
    % Medimos la distancia para ver si hay que aumentar la resolucion
    Distancia = sqrt((1 + DistanciaX - cos(Thetha))^2 + sin(Thetha)^2);
    if (Distancia < DistanciaSiguienteReduccion)
        DistanciaSiguienteReduccion = DistanciaSiguienteReduccion / CofReductorDistancia;
        NumIterThetha = NumIterThetha * CofIncremNumIterThetha;
        NumIterRho = NumIterRho * CofIncremNumIterRho;
        ResolRho = (2*pi) / NumIterRho;
        ResolThetha = pi / NumIterThetha;
    end
    
    % La unidad de superficie indica por cuanta carga computa cada trozito
    % de carga que estamos calculando y su suma debería valer 4*pi
    UnidadSuperficie = abs(sin(Thetha)) * ResolRho * ResolThetha;
    
    Rho = ResolRho/2;   
    while(Rho <= (2*pi))
    
        Pos = [cos(Thetha), abs(sin(Thetha))*sin(Rho), abs(sin(Thetha))*cos(Rho)];
        VectorOrigienDestino = PosicionEvaluada - Pos;
        
        [E, B] = CampoRadiadoCargaAcelerada2(VectorOrigienDestino, Velocidad, Aceleracion);
        E_rad_tot = E_rad_tot + E * UnidadSuperficie;
        B_rad_tot = B_rad_tot + B * UnidadSuperficie;
        
        if( norm(E) > norm(B)*1.0001)
            disp('mierda')
        end
        
        Rho = Rho + ResolRho; 
        num_iteraciones_hechas = num_iteraciones_hechas + 1;
        superficie_acumulada = superficie_acumulada + UnidadSuperficie;
    end
    
   Thetha = Thetha - ResolThetha;
end

% Mostrar resultados
num_iteraciones_hechas = num_iteraciones_hechas
superficie_acumulada = superficie_acumulada

[E_rad_centro, B_rad_centro] = CampoRadiadoCargaAcelerada2(PosicionEvaluada, Velocidad, Aceleracion);

E_rad_centro = E_rad_centro * superficie_acumulada
B_rad_centro = B_rad_centro * superficie_acumulada

E_rad_tot = E_rad_tot
B_rad_tot = B_rad_tot 

% Ratio_E_rad_tot_div_E_rad_centro = norm(E_rad_tot) / norm(E_rad_centro)
% Ratio_B_rad_tot_div_B_rad_centro = norm(B_rad_tot) / norm(B_rad_centro)

S = cross(E_rad_tot, B_rad_tot)

