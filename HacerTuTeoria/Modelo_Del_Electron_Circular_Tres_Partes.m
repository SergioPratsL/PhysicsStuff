%%% 31.03.2023
%syms A q1 q2 m1 m2;

% m1 = 0.5 + 0.25*A;
% m2 = 0.25 - 0.125*A;
% q1 = 0.5 + 0.5*A;
% q2 = 0.25 - 0.25*A;
% 
% 
% f = 2*q1*q2/m1 - (q2^2+q1*q2)/m2
% 
% %f = - ((A/2 + 1/2)*(A/4 - 1/4) - (A/4 - 1/4)^2)/(A/8 - 1/4) - ((A + 1)*(A/4 - 1/4))/(A/4 + 1/2)
% 
% A = 1
% valor = subs(f)
% 
% eqn = f == 0;
% S = solve(eqn);

% Lo siguiente es un intento desesperado

% m1 = 0.5*A;
% m2 = 0.5 - 0.25*A;
% q1 = A;
% q2 = 0.5 - 0.5*A;
% 
% 
% f = 2*q1*q2/m1 - (q2^2+q1*q2)/m2
% 
% A = 1
% valor = subs(f)
% 
% eqn = f == 0;
% S = solve(eqn)


% 07.04.2023. Probar con 3 variables diferentes
denom_q1_ini = 1;
num_q1_ini = -5;
denom_m1_ini = 1;
num_m1_ini = -5;


Cofs = [1, 1, 0, 0; 0, 0, 1, 1; 1, -1, 0, 0; 0, 0, 1, -1];
%A = 1;
A = -1;
%A = 1/2;
%A = -1/2;
%A = 2;
%A = -2
cont = 0;
resultados = 0;

denom_q1 = denom_q1_ini;
num_q1 = num_q1_ini;
denom_m1 = denom_m1_ini;
num_m1 = num_m1_ini;

% Ataque por fuerza bruta
while denom_q1 < 10
    while num_q1 < 6
        while denom_m1 < 10
            while num_m1 < 6
                m1 = num_m1 / denom_m1;
                q1 = num_q1 / denom_q1;
                
                Valores = [1-m1, 1-q1, A/2-m1, A-q1];
                X = linsolve(Cofs, Valores');
                
                m2 = X(1);
                m3 = X(2);
                q2 = X(3);
                q3 = X(4);
                
                dif1 = (q1*q2+q1*q3)*m2 - (q1*q2+q2*q3)*m1;
                dif2 = (q1*q2+q1*q3)*m3 - (q1*q3+q2*q3)*m1;
                
                if abs(dif1) < 0.0000001 && abs(dif2) < 0.0000001 
                    valid_sol = [m1, m2, m3, q1, q2, q3];
                    resultados = resultados + 1;
                    disp(valid_sol)
                end
                
                num_m1 = num_m1 + 1;
                if num_m1 == 0
                    num_m1 = 1;
                end
                cont = cont + 1;
            end
            num_m1 = num_m1_ini;
            denom_m1 = denom_m1 + 1;
        end
        denom_m1 = denom_m1_ini;
        num_q1 = num_q1 + 1;
        if num_q1 == 0
            num_q1 = 1;
        end        
    end
    num_q1 = num_q1_ini;
    denom_q1 = denom_q1 + 1;    
end

disp('done')
cont
resultados

