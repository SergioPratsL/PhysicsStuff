% La idea es probar esto con Schrödinger para ver que va bien y entonces
% intentarlo con Dirac.
% https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator
clear all; clc;

% --- Physical Constants (SI Units) ---
m_e = 9.10938356e-31;     % Electron mass (kg)
c = 299792458;          % Speed of light (m/s)
h_bar = 1.0545718e-34;   % Reduced Planck constant (J.s)
q = 1.602176 * 10^-19;

mc2 = m_e * c^2; 
hc = h_bar * c;

% --- Simulation Parameters ---

%k= -1;  % NO USAR NEGATIVOS
%k = 1;          
%k=2;
%k=3;
%k=4;             
%k=8;   % Comienza a ir mal          
%k = 10;
k=12;  % Revienta y con relleno también
%k=16;  % Con relleno = 4 acaba reventando igualmente.

% En la wiki la energía es k/2 y así lo haré aquí.
%w_wiki = sqrt(k/m_e);
w_wiki = sqrt(abs(k)/m_e);
nivel_base = h_bar * w_wiki / 2;

E = nivel_base + mc2;  
%E = 0.5 * nivel_base + mc2;  
%E = 0.99 * nivel_base + mc2;  

Sigma_Gaussiana = m_e*w_wiki/(2*h_bar)

Sx = MatrizSpin_4_4(1);
Sy = MatrizSpin_4_4(2);
Sz = MatrizSpin_4_4(3);
[at, ax, ay, az] = MatricesAlfa();

% Con estos valores sólo puedo calcular n's pares, pero ya está bien.
x0 = 0;
% Solucion que tiene phi_A no nulo, orientada en Z
%f0 = [1, 0, 0, 0];

% Solucion de spin en X!!!
%f0 = [1, 1, 0, 0] / sqrt(2);

% Ángulo de 45%
f0 = [cos(pi/8), sin(pi/8), 0, 0];

% Numerical settings
%num_iterations = 10000;
% Más iteraciones que ha de llegar mucho más lejos!
num_iterations = 100000;
dx = 1e-14;             % Step size (m)
% Pasos más grandes
%dx = 2*1e-14;             % Step size (m)
    
% Initialize variables for this run
x = zeros(num_iterations, 1);
f = zeros(num_iterations, 4);
    
x(1) = x0;
f(1, 1:4) = f0;  

f_prim = zeros(1,4);
max_change_per_step = 0.01;

% Pre-calculate constants outside the loop
inv_hc = -1i/hc;
x = (0:num_iterations-1) * dx;

disp(['Evaluando la energia ', num2str(E), ' para soluciones pares']);

%factor_relleno = 4;
factor_relleno = 1;
num_iterations = num_iterations * factor_relleno;
dx = dx / factor_relleno;
n_relleno = 1;
j = 2;

dxSquare = dx^2;
f_old = f0;

for i = 2:num_iterations-1

    %x(i) = x(i-1) + dx;
    pot = 0.5 * k *dxSquare*(i-1)^2;
    
    %f_now = f(i-1, :);
    
    f_prim(1) = inv_hc * (E + mc2 - pot)*f_old(4);
    f_prim(2) = inv_hc * (E + mc2 - pot)*f_old(3);
    f_prim(3) = inv_hc * (E - mc2 - pot)*f_old(2);
    f_prim(4) = inv_hc * (E - mc2 - pot)*f_old(1);
    
    % En este punto ya no es viejo, pero prefieron no usar otra variable
    f_old = f_old + f_prim*dx;
    
    if n_relleno == factor_relleno
        n_relleno = 1;
        f(j,1:4) = f_old; 
        j = j + 1;
    else
        n_relleno = n_relleno + 1;
    end
end

j = j - 1;  % Quitar el último incremento que no queremos.

rho = sum(abs(f).^2, 2);
    
disp(['Distancia al centro: ', num2str(x(j))]);

valor_exponente_gaussiana_final = - m_e*w_wiki*x(j)^2 / (2*h_bar);

% 1. Compute Sx * f' (Apply the operator to all spinors)
% 2. Multiply by f (dot product for each row)
% 3. Normalize by rho
SpinZ = sum((f * Sz.') .* conj(f), 2) ./ rho;
SpinY = sum((f * Sy.') .* conj(f), 2) ./ rho;
SpinX = sum((f * Sz.') .* conj(f), 2) ./ rho;

%Jx = sum((f * ax.') .* conj(f), 2) ./ rho;
%Jy = sum((f * ay.') .* conj(f), 2) ./ rho;
%Jz = sum((f * ay.') .* conj(f), 2) ./ rho;

%rho_large = sum(abs(f(:, 1:2)).^2, 2);
%rho_small = sum(abs(f(:, 3:4)).^2, 2);


iterations = 1:(num_iterations/factor_relleno);

figure;

% Top Plot: Rho
subplot(3,1,1); 
plot(iterations, rho, 'LineWidth', 1.5);
ylabel('\rho (Magnitude)');
grid on;
title('Bispinor Evolution');

% Middle Plot: SpinZ
subplot(3,1,2);
%plot(iterations, rho_large, 'r', 'LineWidth', 1.5);
%ylabel('\langle rho_large \rangle');
plot(iterations, SpinZ, 'r', 'LineWidth', 1.5);
ylabel('\langle S_z \rangle');
%plot(iterations, Jx, 'r', 'LineWidth', 1.5);
%ylabel('\langle Jx \rangle');
xlabel('Iteration Index');
grid on;

% Bottom Plot: SpinY
subplot(3,1,3);
%plot(iterations, rho_small, 'r', 'LineWidth', 1.5);
%ylabel('\langle rho_small \rangle');
%plot(iterations, SpinY, 'r', 'LineWidth', 1.5);
%ylabel('\langle S_y \rangle');
plot(iterations, SpinX, 'r', 'LineWidth', 1.5);
ylabel('\langle S_x \rangle');
%plot(iterations, Jy, 'r', 'LineWidth', 1.5);
%ylabel('\langle Jy \rangle');
xlabel('Iteration Index');
grid on;


% El primer test, con k=1 da esto:
% f(25000, :) = 0.7537 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 - 0.0003i
% f(50000, :) = 0.3226 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 - 0.0003i
% f(75000, :) = 0.0785 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 - 0.0001i
% f(99000, :) = 0.0121 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 - 0.0000i


