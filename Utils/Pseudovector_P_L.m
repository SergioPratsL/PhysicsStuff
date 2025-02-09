function W = Pseudovector_P_L( M, P )
% J es el tensor relativista de momento angular
% P es el 4 momento

W = [0,0,0,0];

% Implementación de Levi Civita con preciosos bucles for :D
for i=1:4
    
    for j=1:4
        if j == i
            continue
        end
        
        for k=1:4
            if k == j || k == i
                continue
            end
            
            for m=1:4            
                if m == k || m == j || m == i
                    continue
                end     
                
%                 if ((j==1 && k == 2) || (j==2 && k == 1))
%                     a = 1;
%                 end         
                
                cof = LeviCivita([i,j,k,m]);
                
                W(i) = W(i) + cof * M(j,k) * P(m);
            end
        end
    end
end

% Finalmente aplico el factor 1/2
W = W / 2;
