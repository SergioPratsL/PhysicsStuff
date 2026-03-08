% Como el QHO_Dirac_test pero con E cte bueno E depende del signo de X, es
% como su en el 0 hubiera un plano con densidad superficial de corriente.

clear all; clc;

% --- Physical Constants (SI Units) ---
m = 9.10938356e-31;     % Electron mass (kg)
c = 299792458;          % Speed of light (m/s)
hbar = 1.0545718e-34;   % Reduced Planck constant (J.s)
q = 1.602176 * 10^-19;

% --- Simulation Parameters ---
% k incluya q de algumna manera...
%k = 1000;               % QHO Spring constant (N/m) - Example value
% Más caña para que llegue antes a la energía negativa!
%k = 10000;
% menos caña
k = 0.001;

E_mc2 = m * c^2; 

% Si E=mc2 y k=0, g no levanta cabeza y f está plano, lógico.
%k=0;

% En la wiki la energía es k/2, yo uso k, por ello el su k es el doble que
% el mio
w_wiki = sqrt(2*k/m);
nivel_base_old = hbar * w_wiki / 2;

nivel_base = sqrt(E_mc2^2 + 2*k*c*hbar) - E_mc2;

% --- Energy Scan Settings ---
energy_iterations = 100;
%energy_iterations = 1;

% Defined Interval for Energy Scan (E_min to E_max)
% Scanning around the rest mass energy + some oscillator levels
% Adjust these factors to widen/narrow the search
%E_min = E_mc2 + 0.1 * nivel_base;  

%E_min = E_mc2;  
E_min = E_mc2 + 0.45 * nivel_base;  
%E_max = E_mc2 + 2 * nivel_base; 
E_max = E_mc2 + 0.6 * nivel_base; 
energies = linspace(E_min, E_max, energy_iterations);

% Initial conditions at x = 0
x0 = 0;
f0 = 1;                 % f(0) = 1
f_prime_0 = 0;          % f'(0)

% Numerical settings
num_iterations = 10000;
dx = 1e-14;             % Step size (m)
% Pasos más grandes
%dx = 2*1e-14;             % Step size (m)
tail_size = num_iterations / 10;       % Last 10% steps to measure probability

% --- Initialization ---
scan_prob_results = zeros(energy_iterations, 1);
min_wave_prob = Inf;
best_E_index = -1;

disp(['Starting Energy Sweep over ', num2str(energy_iterations), ' values...']);
filename = 'E_cte_test_results.txt';
fileID = fopen(filename, 'w');
fprintf(fileID, 'Iter\t Energy(J)\t Mean_Tail_Prob\n');

% --- Outer Loop: Energy Sweep ---
for j = 1:energy_iterations
    E = energies(j);
    
    % Initialize variables for this run
    x = zeros(num_iterations, 1);
    f = zeros(num_iterations, 1);
    g = zeros(num_iterations, 1);
    
    x(1) = x0;
    f(1) = f0;
    
    % Calculate g(0) for current E
    g(1) = (hbar * c * f_prime_0) / (E + m * c^2);
    
    current_tail_prob_sum = 0;
    start_tail_idx = num_iterations - tail_size;
    
    % --- Inner Loop: Spatial Integration ---
    for i = 1:num_iterations-1
        xi = x(i);
        
        E_minus = (E - m*c^2 - k*xi);
        E_plus = (E + m*c^2 - k*xi);
        
        % Derivatives
        dg_dx = (E_minus * f(i)) / (hbar * c);
        df_dx = (E_plus * g(i)) / (hbar * c);
        
        if (i == 100)
            sdasd = 1;
        end
        
        % Update
        x(i+1) = xi + dx;
        f(i+1) = f(i) + dx * df_dx;
        g(i+1) = g(i) + dx * dg_dx;
        
        % Accumulate Probability in Tail
        if i >= start_tail_idx
            current_tail_prob_sum = current_tail_prob_sum + (f(i+1)^2 + g(i+1)^2);
        end
    end
    
    % Store Scan Result
    avg_tail_prob = current_tail_prob_sum / tail_size;
    scan_prob_results(j) = avg_tail_prob;
    fprintf(fileID, '%d\t %.6e\t %.6e\n', j, E, avg_tail_prob);
    
    if current_tail_prob_sum < min_wave_prob
        min_wave_prob = current_tail_prob_sum;
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

