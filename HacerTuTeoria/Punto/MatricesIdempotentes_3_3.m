
vector_inicial = [0,0,0,0,0,0,0,0,0];

%enoontrados = ValoresAProbar(vector_inicial, 1, 0)

% Aparecieron 384 matrices idempotentes, pero muchas serán  fases, ahora
% habría que ver de todas estas matrices cuales tienen

% Mejor mirar esto:
% https://math.stackexchange.com/questions/1728712/higher-dimensional-pauli-matrices

w = cos(pi/3) + 1i * sin(pi/3);

M1 = [0,1,0; 0,0,1; 1,0,0];
M2 = [0,0,w^2; w,0,0; 0,1,0];

M1 * M1';

M2 * M2';

M1*M2' + M2*M1';

% Nope.

% Cosas que me dijo Gemini de Google:

% Pauli matrices are essentially 4-dimensional split-complex numbers. Split-complex numbers can be only 2n-dimensional.

J1 = [0,0,0; 0,0, i; 0, i, 0];

avd = J1 * J1;

% A truly linear Hamiltonian in the same way as Dirac's isn't possible for a scalar (spin-0) relativistic particle. Here's why and what the alternatives are:

 % Spin-0 particles transform differently under Lorentz transformations (they are scalars), 
 % and this fundamentally changes the mathematical requirements.  You can't construct 
 % a similar set of matrices that would allow for a linear Hamiltonian for a scalar particle in the same manner.
 % While you can rewrite it as a first-order system, it involves introducing additional components to the wavefunction, effectively making it not a single-component scalar anymore.
 
%The idea is to get a Hamiltonian based only in first derivatives for spinless particles, the same way as the Dirac Hamiltonian is based on first derivarives for particles with spin

% In Summary:

% You can't create a truly analogous, single-component, first-order Hamiltonian for spin-0 particles like the Dirac Hamiltonian. 
% The fundamental reason lies in the different Lorentz transformation properties of spin-0 and spin-1/2 particles. 
% The Klein-Gordon equation, while second-order, is the correct and fundamental equation for relativistic spin-0 particles. 
% While first-order reformulations exist, they involve additional components and aren't a direct parallel to the Dirac equation's single-component linearization.

% La reformulación que se ofrece es crear una variable que es de primer
% grado en tiempo pero no en espacio, es muy cutre, pero parece que no se
% puede...

function encontrados = ValoresAProbar(vector_actual, posicion_actual, encontrados_ini)

    valores_a_probar = [0, 1, -1, 1i, -1i];

    num_valores = 5;
    n = 1;

    encontrados = encontrados_ini;

    while n <= num_valores

        vector_actual(posicion_actual) = valores_a_probar(n);

        if posicion_actual < 9
            encontrados = ValoresAProbar(vector_actual, (posicion_actual + 1), encontrados);
        else
            encontrados = encontrados + EvaluaMatriz(vector_actual);
        end

        n = n + 1;
    end    

end


function encontrado = EvaluaMatriz(vector_actual)
   
    encontrado = 0;
    
    matriz = [vector_actual(1), vector_actual(2), vector_actual(3); vector_actual(4), vector_actual(5), vector_actual(6); vector_actual(7), vector_actual(8), vector_actual(9)];
    
    producto = matriz * matriz';
    
    if( producto(1,2) ~= 0 || producto(1,3) ~= 0 || producto(2,1) ~= 0 || producto(2,3) ~= 0 || producto(3,1) ~= 0 || producto(3,2) ~= 0)
        return
    end
    
    if ( producto(1,1) == 0 || producto(2,2) ~= producto(1,1) || producto(3,3) ~= producto(1,1))
        return
    end
    
    encontrado = 1;
    
    matriz_idempotente = matriz
    
end

