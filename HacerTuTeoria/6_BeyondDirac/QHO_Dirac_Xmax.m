% QHO_Dirac_test.m
% Simulation of 1D Dirac equation with Harmonic Oscillator potential
% En este caso el potencial no crece hasta el infinito sino que se detiene
% al llegar a xmax, que es donde comienza la simulacion
clear all; clc;

% --- Physical Constants (SI Units) ---
m = 9.10938356e-31;     % Electron mass (kg)
c = 299792458;          % Speed of light (m/s)
hbar = 1.0545718e-34;   % Reduced Planck constant (J.s)
q = 1.602176 * 10^-19;

% --- Simulation Parameters ---
% Adjust these values as needed
% k = 1000;               % QHO Spring constant (N/m) - Example value
% Más caña para que llegue antes a la energía negativa!
k = 10000;

E_mc2 = m * c^2; 

% En la wiki la energía es k/2, yo uso k, por ello el su k es el doble que
% el mio
w_wiki = sqrt(2*k/m);
nivel_base = hbar * w_wiki / 2;

% --- Energy Scan Settings ---
energy_iterations = 100;
%energy_iterations = 1;

% Defined Interval for Energy Scan (E_min to E_max)
% Scanning around the rest mass energy + some oscillator levels
% Adjust these factors to widen/narrow the search
%E_min = E_mc2 + 0.1 * nivel_base;  

E_min = E_mc2 - 0.5* nivel_base;  
%E_max = E_mc2 + 1 * nivel_base; 
E_max = E_mc2 + 0.5 * nivel_base; 
energies = linspace(E_min, E_max, energy_iterations);

% Initial conditions at x = 0
x0 = 0;
f0 = 1;                 % f(0) = 1
f_prime_0 = 0;          % f'(0)

% Numerical settings
%num_iterations = 10000;
% Más iteraciones que ha de llegar mucho más lejos!
num_iterations = 100000;
dx = 2*1e-14;             % Step size (m)
% Pasos más grandes
%dx = 2*1e-14;             % Step size (m)

x_max = (num_interations-1) * dx;

% --- Initialization ---
scan_prob_results = zeros(energy_iterations, 1);
%min_wave_prob = Inf;
min_error = Inf;
best_E_index = -1;

disp(['Starting Energy Sweep over ', num2str(energy_iterations), ' values...']);
filename = 'QHO_Xmax_test_results.txt';
fileID = fopen(filename, 'w');
fprintf(fileID, 'Iter\t Energy(J)\t Mean_Tail_Prob\n');

% --- Outer Loop: Energy Sweep ---
for j = 1:energy_iterations
    E = energies(j);
    
    k1 = E - E_mc2 - k * x_max^2;
    k2 = E + E_mc2 - k * x_max^2;
    w = sqrt(k1/k2);
    
    % Initialize variables for this run
    x = zeros(num_iterations, 1);
    f = zeros(num_iterations, 1);
    g = zeros(num_iterations, 1);

    x(1) = x_max;
    
    if (w < 0)
        f(1) = 1;
        g(1) = w;
    else
        f(1) = 1/w;
        g(1) = 1;
    end
    
    % No hace falta almacenar todas las f_prime y g_prime.
    df_dx = f(1) * k1*k2;
    dg_dx = g(1) * k1*k2;    
   
    
    % --- Inner Loop: Spatial Integration ---
    for i = 1:num_iterations-1       
        xi = x(i);
        E_minus = (E - m*c^2 - k*xi^2);
        E_plus = (E + m*c^2 - k*xi^2);
        
        % Derivatives
        dg_dx = (E_minus * f(i)) / (hbar * c);
        df_dx = (E_plus * g(i)) / (hbar * c); 
        
        if (i == num_iterations-1)
            asad = 1;
        end
        
        % Recordemos que vamos hacia atrás, por eso las derivadas restan
        x(i+1) = xi - dx;
        f(i+1) = f(i) - dx * df_dx;
        g(i+1) = g(i) - dx * dg_dx;
    end
    
    % Valor normalizado
    integral_dens_onda = (norm(f)^2 + norm(g)^2) / num_iterations;
    
    factor_ori = sqrt((E-E_mc2)/(E+E_mc2));
    
    % Estas fórmulas me dan mucho miedo.
    error_1 = ((E-E_mc2)*f(num_iterations))^2 + (dg_dx)^2;
    error_2 = ((E+E_mc2)*g(num_iterations))^2 + (df_dx)^2;
    
    error = min(error_1, error_2);
    
    fprintf(fileID, '%d\t %.6e\t %.6e\n', j, E, error);
    
    if error < min_error
        min_error = error;
        best_E_index = j;
    end
end

fclose(fileID);

% Results Summary
if (best_E_index > -1)
    best_E_val = energies(best_E_index);
    best_avg_prob = min_wave_prob / tail_size;
else
    % The best_E_index was a Nan
    best_E_val = -1;
    best_avg_prob = 1;
end

disp('--------------------------------------------------');
disp('Sweep Complete.');
disp(['Iteration with least wave probability: ', num2str(best_E_index)]);
disp(['Corresponding Energy: ', num2str(best_E_val), ' J']);
disp(['Average Tail Probability (sum/1000): ', num2str(best_avg_prob)]);
disp(['Results summary saved to ', filename]);
disp('--------------------------------------------------');

% --- Simple Plot ---
figure;
subplot(2,1,1);
plot(x, f, 'b', 'LineWidth', 1.5);
ylabel('f(x)'); title(['QHO Dirac Solution (E = ', num2str(E), ' J)']);
grid on;

subplot(2,1,2);
plot(x, g, 'r', 'LineWidth', 1.5);
xlabel('x (m)'); ylabel('g(x)');
grid on;
