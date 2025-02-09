
% El objetivo de estos scripts es calcular la proyección de estados de
% momento angula en X y en Y.

NumIterThetha = 500;
NumIterPhi = 500;

ResolPhi = (2*pi) / NumIterPhi;
ResolThetha = pi / NumIterThetha;


% Empezamos por Theta cercano a 180º
ThethaInicial = pi - ResolThetha/2;
Thetha = ThethaInicial;

SuperficieAcumulada = 0;
NumIteracionesHechas = 0;
ProdTot = 0;

while Thetha > 0   
    % La unidad de superficie indica por cuanta carga computa cada trozito
    % de carga que estamos calculando y su suma debería valer 4*pi
    UnidadSuperficie = abs(sin(Thetha)) * ResolPhi * ResolThetha;
    
    Phi = ResolPhi/2;   
    while(Phi <= (2*pi))
        % Pos tiene la posición que estamos evaluando, para los estados
        % cuyo eje es la X es esencial.
        Pos = [cos(Thetha), abs(sin(Thetha))*sin(Phi), abs(sin(Thetha))*cos(Phi)];
        
        Prod = cos(Thetha)^2;
        
        ProdTot = ProdTot + Prod * UnidadSuperficie;
        
        Phi = Phi + ResolPhi; 
        NumIteracionesHechas = NumIteracionesHechas + 1;
        SuperficieAcumulada = SuperficieAcumulada + UnidadSuperficie;
    end
    
   Thetha = Thetha - ResolThetha;
end

% Da 12.5664 que es 4pi :)
%SuperficieAcumulada = SuperficieAcumulada

FactorAmplitud = 1/sqrt(2*pi) * sqrt(6)/2;

% Da 1 :)
ProdTot = ProdTot * FactorAmplitud^2
