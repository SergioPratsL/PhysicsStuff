% La idea es probar esto con Schrödinger para ver que va bien y entonces
% intentarlo con Dirac.
% https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator
clear all; clc;

% --- Physical Constants (SI Units) ---
m_e = 9.10938356e-31;     % Electron mass (kg)
%c = 299792458;          % Speed of light (m/s)
h_bar = 1.0545718e-34;   % Reduced Planck constant (J.s)
%q = 1.602176 * 10^-19;

% --- Simulation Parameters ---

%k=10000
%k = 1000;
%k = 100;

%k = -10000;
%k = -1000;
%k = -100;
%k= -10;
%k= -1;
%k = 1;          % El mejor y va mal si E=2*nivel_base!
%k=2;            % Parece una gaussiana y quiere E correcto
%k=4;             % Es poca cosa, pero empieza a ir mal   
%k=6;             % El desastre crece
k=10;
k=8;

% En la wiki la energía es k/2 y así lo haré aquí.
%w_wiki = sqrt(k/m_e);
w_wiki = sqrt(abs(k)/m_e);
nivel_base = h_bar * w_wiki / 2;

E = nivel_base;  
%E = 2 *nivel_base;  
%E = 10 * nivel_base;

Sigma_Gaussiana = m_e*w_wiki/(2*h_bar)

% Con estos valores sólo puedo calcular n's pares, pero ya está bien.
x0 = 0;
f0 = 1;                 % f(0) = 1
f_prime_0 = 0;          % f'(0)

% Numerical settings
%num_iterations = 10000;
% Más iteraciones que ha de llegar mucho más lejos!
num_iterations = 100000;
dx = 1e-14;             % Step size (m)
% Pasos más grandes
%dx = 2*1e-14;             % Step size (m)
    
% Initialize variables for this run
x = zeros(num_iterations, 1);
f = zeros(num_iterations, 1);
    
x(1) = x0;
f(1) = f0;
f_prim = 0;   

max_change_per_step = 0.01;

disp(['Evaluando la energia ', num2str(E), ' para soluciones pares']);
   
for i = 2:num_iterations-1

    x(i) = x(i-1) + dx;
        
    f(i) = f(i-1) + f_prim*dx; 
    
    if (i == 95000)
        asda = 10;
    end
    
    if (f(i) ~= 0)
        f_prim2 = -(E - 0.5*k*x(i-1)^2)*f(i-1)*(2*m_e/h_bar^2);
    else
        f_prim2 = 0;
    end
    
    f_prim = f_prim + f_prim2 * dx;
    
    if (f_prim == 0 && f(i) == 0)
        disp("Tanto f como f_prim son 0 ¡Estamos perdidos!")
    end
   
    %if (abs(f_prim)*dx > max_change_per_step * abs(f(i)) && abs(f_prim2)*dx > max_change_per_step * abs(f_prim))
    if (abs(f_prim2)*dx > max_change_per_step * abs(f_prim) && i > 500)
        disp(['Número de iteraciones ', num2str(i)]);
        disp(['Distancia al centro: ', num2str(x(i))]);
        break;
    end
    
end
    
disp(['Distancia al centro: ', num2str(x(i))]);

valor_exponente_gaussiana_final = - m_e*w_wiki*x(i)^2 / (2*h_bar)

% --- Simple Plot ---
figure;
plot(x, f, 'b', 'LineWidth', 1.5);
ylabel('f(x)'); title(['QHO Schro Solution (E = ', num2str(E), ' J)']);
grid on;
