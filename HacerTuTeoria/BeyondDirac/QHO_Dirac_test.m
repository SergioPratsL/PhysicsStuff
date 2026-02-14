% QHO_Dirac_test.m
% Simulation of 1D Dirac equation with Harmonic Oscillator potential
clear all; clc;

% --- Physical Constants (SI Units) ---
m = 9.10938356e-31;     % Electron mass (kg)
c = 299792458;          % Speed of light (m/s)
hbar = 1.0545718e-34;   % Reduced Planck constant (J.s)

% --- Simulation Parameters ---
% Adjust these values as needed
k = 1000;               % QHO Spring constant (N/m) - Example value

% Si E=mc2 y k=0, g no levanta cabeza y f está plano, lógico.
%k=0;

% En la wiki la energía es k/2, yo uso k, por ello el su k es el doble que
% el mio
w_wiki = sqrt(2*k/m);
nivel_base = hbar * w_wiki / 2;

% Energy E must be > m*c^2
% Los dos caen a valores negativos muy grandes.
%E = m * c^2 + 1.0e-18;  % Total Energy (Rest mass + kinetic/potential)

% nivel base del oscilador cuántico clásico
% Los dos oscilan, pero la oscilación va a más!
%E = m * hbar * nivel_base;

% Se van de madre hacia lo negativo, f parece frenar...
%E = m * c^2;

% Oscilaciones amplificadas.
E = m * hbar * 2 * nivel_base;


% Initial conditions at x = 0
x0 = 0;
f0 = 1;                 % f(0) = 1
f_prime_0 = 0;          % f'(0) - Determining g(0)
                        % For even parity solutions (Ground state), f'(0) is typically 0.

% Numerical settings
num_iterations = 10000;
dx = 1e-14;             % Step size (m). Keep small for atomic scales!

% --- Initialization ---
x = zeros(num_iterations, 1);
f = zeros(num_iterations, 1);
g = zeros(num_iterations, 1);

x(1) = x0;
f(1) = f0;

% Initial calculation of g(0) from f'(0)
% Eq 2: (E + m*c^2 - k*x^2)g = hbar*c*df/dx
% At x=0: g(0) = (hbar*c * f'(0)) / (E + m*c^2 - k*x(0)^2)
g(1) = (hbar * c * f_prime_0) / (E + m * c^2 - k * x(1)^2);

% --- Iteration Loop ---
disp('Starting simulation...');
filename = 'QHO_test_results.txt';
fileID = fopen(filename, 'w');
fprintf(fileID, 'Index\t x(m)\t f(x)\t g(x)\n');

% Print initial state
fprintf(fileID, '%d\t %.6e\t %.6e\t %.6e\n', 1, x(1), f(1), g(1));

for i = 1:num_iterations-1
    % Current x
    xi = x(i);
    
    % Calculate derivatives at current step
    % Eq 1: (E - m*c^2 - k*x^2)f = hbar*c * dg/dx
    % dg/dx = [(E - m*c^2 - k*x^2) * f] / (hbar*c)
    dg_dx = ((E - m*c^2 - k*xi^2) * f(i)) / (hbar * c);
    
    % Eq 2: (E + m*c^2 - k*x^2)g = hbar*c * df/dx
    % df/dx = [(E + m*c^2 - k*x^2) * g] / (hbar*c)
    df_dx = ((E + m*c^2 - k*xi^2) * g(i)) / (hbar * c);
    
    % Update next values (Euler method)
    x(i+1) = xi + dx;
    f(i+1) = f(i) + dx * df_dx;
    g(i+1) = g(i) + dx * dg_dx;
    
    % Write to file (optional: write every N steps to save space if needed)
    fprintf(fileID, '%d\t %.6e\t %.6e\t %.6e\n', i+1, x(i+1), f(i+1), g(i+1));
    
    % Check stop condition for turning point
    % We want to see if k*x^2 > E + m*c^2
    % This is the "hard" turning point for the g component.
    % The classical allowable region is typically E > m*c^2 + 0.5*k*x^2 (approx).
    limit_val = E + m*c^2;
    potential_term = k * x(i+1)^2;
    
    if i == num_iterations-1
        fprintf('Reached max iterations. x_final = %.3e\n', x(i+1));
        if potential_term <= limit_val
            warning('Did not reach the condition k*x^2 > E + m*c^2. Increase num_iterations or dx.');
        end
    end
end

fclose(fileID);
disp(['Simulation complete. Results saved to ', filename]);

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
