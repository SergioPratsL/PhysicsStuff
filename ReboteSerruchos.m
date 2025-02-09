function p_x_fin = RebotaSerruchos( pos_ini, Rho_ini, Psi_ini, H, L, tan_alfa, cos_alfa , sen_alfa )
% Pos contiene [x,y] (z es irrelevante).
% Psi_ini necesario porque Psi cambia al rebotar con el muro.


p_x = 0;
Rebota = 1;
pos = pos_ini;
Rho = Rho_ini;
Psi = Psi_ini;

% Pendiente expresada como unidades Y avanzadas por cada unidad X (ignorando variacion en Z).    
if mod(Psi, pi/2) >= pi / 2  - 0.0001 || mod(Psi, pi/2) <= 0.0001
    if ( Rho < 0 && Psi <= pi / 2 )  || ( Rho > 0 && Psi > pi/2 )
       pendiente = -10000;          % Aproximacion :)
   else
        pendiente = 10000;          % Aproximacion :)
    end
else
    pendiente = tan(Rho) / cos(Psi);
end           


if pos(2) > H && Rho >= 0
    Rebota = 0;
    p_x_fin = 0;
end


% Vamos a ver en que punto llega a estar a altura H
if pos(2) > H
    dif_x = (pos(2) - H) / pendiente;
% Actualizamos la posicion
    pos(2) = H;
    pos(1) = mod( (pos(1) + dif_x), L );
    
end

cont_aux = 1;
while Rebota == 1 && cont_aux < 100

% Ahora no tengo mas que entrar en las casuisticas:

% En este caso siempre rebotara en la rampa
    if pendiente <= 0 && Psi <= pi /2
    
% Ecuacion a resolver:
    % pos(2) + (x - pos(1)) * pdte = tan_alfa * x
        x = ( pos(2) - pos(1) * pendiente ) / ( tan_alfa - pendiente);
        pos(1) = x;
        pos(2) = x * tan_alfa;
% Evaluar el impacto en muro    
        [pendiente, dif_p, cambiaPsi] = MomentoRebote(pendiente, 'rampa', tan_alfa, Psi);
        if cambiaPsi == 1
            Psi = - Psi + pi;
        end                              
        
        p_x = p_x + dif_p;
        
% Determinar si rebota o escapa        
    elseif pendiente > 0 && Psi < pi/2
        
        if tan_alfa <= pendiente 
% Escapa seguro
            Rebota = 0;
        else
% Solo deberia pasar para Rebotas en el muro por lo que pos(1) = 0 pero lo
% incluyo por si acaso el matlab deja restos del tipo e-15 o asi...
            if (tan_alfa - pendiente) * (L - pos(1)) > (H - pos(2))
% pos(2) + pendiente * x = x * tan_alfa                
                x = pos(2) / (tan_alfa - pendiente);
                pos(2) = x * tan_alfa; 
                pos(1) = x;
                [pendiente, dif_p, cambiaPsi] = MomentoRebote(pendiente, 'rampa', tan_alfa, Psi);
                if cambiaPsi == 1
                    Psi = - Psi + pi;
                end                      
                p_x = p_x + dif_p;       
            else
                Rebota = 0;                       
            end
        end
        
    elseif  Psi >= pi/2 && pendiente < 0  % pendiente < 0 es que sube porque va hacia atras!!
        
        if (pos(2) - pos(1) * pendiente )  > H
            Rebota = 0;
        else
            pos(2) = pos(2) - pos(1) * pendiente;
            pos(1) = 0;
            Psi = - Psi + pi;
            [pendiente, dif_p] = MomentoRebote(pendiente, 'muro_', tan_alfa, Psi);
            p_x = p_x + dif_p;  
        end   
        
    else    % Psi > pi/2 && pendiente >= 0   (realmente baja)
% Si pos(2) - pos(1) * pdte > 0, choca con muro.      

        if  (pos(2) - pos(1) * pendiente) > 0
% En este caso choca contra el muro                
            pos(2) = pos(2) - pos(1) * pendiente;
            pos(1) = 0;
            Psi = - Psi + pi;
            [pendiente, dif_p] = MomentoRebote(pendiente, 'muro_', tan_alfa, Psi);                
            p_x = p_x + dif_p;
        else
            x = (pos(2) - pos(1) * pendiente) / (tan_alfa - pendiente);
            pos(2) = x * tan_alfa;
            pos(1) = x;
            [pendiente, dif_p, cambiaPsi] = MomentoRebote(pendiente, 'rampa', tan_alfa, Psi);
            if cambiaPsi == 1
                Psi = - Psi + pi;
            end                
            p_x = p_x + dif_p;                              
        end
        

    end
    
% Salvaguardas    
    if pos(2) > H && Rebota == 1
        pos(2) = H;
    end
    if pos(1) > L && Rebota == 1
        pos(1) = mod( pos(1), L);
    end
        
    cont_aux = cont_aux + 1;
end


p_x_fin = p_x;