
% El objetivo de estos scripts es calcular la proyección de estados de
% momento angular en X y en Y.
% Comparar el estado (2,1,0) en Z contra el (2,1,1) en X
% La amplitud es raiz(2)/2.


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
        
        rXY = sqrt(Pos(1)^2+Pos(2)^2);
        %FactorFase1 = Pos(1)/rXY + 1i*Pos(2)/rXY;
        %Factor1 = rXY;
        % Las dos líneas anteriores se pueden simplificar como:
        Factor1 = Pos(1) + 1i*Pos(2);
        
        rZY = sqrt(Pos(2)^2+Pos(3)^2);
        %FactorFase2 = Pos(3)/rZY + 1i*Pos(2)/rZY;
        %FactorThetha2 = rZY;
        % Las dos líneas anteriores se pueden simplificar como:
        
        % Creo que lo si la X pasa a ser la Z, la Z pasa a ser -X!
        Factor2 = -Pos(3) + 1i*Pos(2);
        
        Prod = Factor1 * Factor2';
        % No conjugar hace que calcules 2_1_1_x_2_1_-1x, que también
        % interesa
        %Prod = Factor1 * Factor2;
        
        ProdTot = ProdTot + Prod * UnidadSuperficie;

        Phi = Phi + ResolPhi; 
        %NumIteracionesHechas = NumIteracionesHechas + 1;
        %SuperficieAcumulada = SuperficieAcumulada + UnidadSuperficie;
    end
    
   Thetha = Thetha - ResolThetha;
end

% Da 12.5664 que es 4pi :)
%SuperficieAcumulada = SuperficieAcumulada

FactorAmplitud1 = 1/sqrt(2*pi) * sqrt(3)/2;
FactorAmplitud2 = FactorAmplitud1;

% Da 1 :)
ProdTot = ProdTot * FactorAmplitud1 * FactorAmplitud2
