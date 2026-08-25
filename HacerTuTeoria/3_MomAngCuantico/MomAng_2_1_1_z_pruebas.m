
% El objetivo de estos scripts es calcular la proyección de estados de
% Comparar el estado (2,1,0) en Z contra el (2,1,1) en Y
% Este caso el resultado es imaginaria y no real.


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
        Pos = [sin(Thetha)*sin(Phi), sin(Thetha)*cos(Phi), cos(Thetha)];
        
        % No tocar factor 1, aquí es (2,1,1)!  (excepto para alternar el
        % signo)
        Factor1 = Pos(1) - 1i*Pos(2);
        
        Factor2 = Pos(1);        
        
        Prod = Factor1 * Factor2';
        
        ProdTot = ProdTot + Prod * UnidadSuperficie;

        Phi = Phi + ResolPhi; 
        NumIteracionesHechas = NumIteracionesHechas + 1;
        SuperficieAcumulada = SuperficieAcumulada + UnidadSuperficie;
    end
    
   Thetha = Thetha - ResolThetha;
end

% Da 12.5664 que es 4pi :)
%SuperficieAcumulada = SuperficieAcumulada

FactorAmplitud1 = 1/sqrt(2*pi) * sqrt(6)/2;
FactorAmplitud2 = 1/sqrt(2*pi) * sqrt(3)/2;

% Da 1 :)
ProdTot = ProdTot * FactorAmplitud1 * FactorAmplitud2
