% Top Plot: Rho
subplot(2,1,1); 
plot(iterations, rho, 'LineWidth', 1.5);
ylabel('\rho (Magnitude)');
grid on;
title('Bispinor Evolution');

% Bottom Plot: SpinZ
subplot(2,1,2);
plot(iterations, SpinZ, 'r', 'LineWidth', 1.5);
ylabel('\langle S_z \rangle');
xlabel('Iteration Index');
grid on;
