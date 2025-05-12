% https://en.wikipedia.org/wiki/Covariant_formulation_of_classical_electromagnetism
% https://en.wikipedia.org/wiki/Electromagnetic_stress%E2%80%93energy_tensor

clear;

E = [11, 13, 17];
B = [3, 5, 7];


abc = 11^2 + 13^2 + 17^2 + 3^2 + 5^2 + 7^2;

u = (norm(E)^2 + norm(B)^2) / 2
p = cross(E, B)

tMaxwell = TensorDeMaxwell(E, B);

F = TensorEM(E, B);


T11 = CalculaTensorEstresEnergia(F, 1, 1)
T12 = CalculaTensorEstresEnergia(F, 1, 2)
T13 = CalculaTensorEstresEnergia(F, 1, 3)
T14 = CalculaTensorEstresEnergia(F, 1, 4)

T21 = CalculaTensorEstresEnergia(F, 2, 1)
T31 = CalculaTensorEstresEnergia(F, 3, 1)
T41 = CalculaTensorEstresEnergia(F, 4, 1)

% T22 = CalculaTensorEstresEnergia(F, 2, 2)
% T23 = CalculaTensorEstresEnergia(F, 2, 3)
% T24 = CalculaTensorEstresEnergia(F, 2, 4)
% T22 = CalculaTensorEstresEnergia(F, 2, 2)
% T23 = CalculaTensorEstresEnergia(F, 2, 3)
% T24 = CalculaTensorEstresEnergia(F, 2, 4)


function val = CalculaTensorEstresEnergia(F, alfa, betta)
    eta = MatrizMinkowski();

    Fc = F';
    %Fc = F;
    
    %val_cte = ContraeMatrices(F, F);  
    val_cte = ContraeMatricesEspecial(F, F);
    
    val = 0;

    iters = 0;
    
    ni = 1;
    while ni <= 4
        gamma = 1;
        cof_eta = eta(alfa, ni);
        cof_eta = cof_eta * eta(betta, betta);
        while gamma <= 4            
            addVal = F(ni,gamma) * F(betta,gamma) * cof_eta;
            if addVal ~= 0
                iters = iters + 1;
            end
            val = val + addVal;
            gamma = gamma + 1;
        end

       ni = ni + 1;
    end
    
    val = val - (1/4) * val_cte * eta(alfa, betta);
    
    iters = iters;
end


function val = ContraeMatricesEspecial(T1, T2)
    val = 0;

    m = 1;
    while m <= 4
        n = 1;
        while n <= 4
            if (m > 1 && n > 1)
                factor = -1;
            else
                factor = 1;
            end
            val = val + factor * T1(m,n) * T2(m,n);
            n = n + 1;
        end

        m = m + 1;
    end
end

