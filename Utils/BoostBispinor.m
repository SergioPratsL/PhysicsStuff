function BispinorBoosteado = BoostBispinor(bispinor, u)
% Aplica boosts a bispinors de Dirac
% 18.10.2026 Esto parece que no es que cambie el observador de la 
% misma onda sino que empuja el boost a velocidad v, esto es absurdo
% y acaba aquí!


if (norm(u) == 0)
    BispinorBoosteado = bispinor;
    return
end

% 18.10.2026. Forma más fácil de arreglar la ñapa.
v = -u;

v_norm = v / norm(v);
rapidez = atanh(norm(v));
rapidez_vect = rapidez * v_norm;

expM = MatrizBoostSpinDirac(rapidez_vect);

%BispinorBoosteado = expM * bispinor';
% Gemini 25.10.2026
BispinorBoosteado = expM * bispinor.';


% Las pruebas parece que dan peor resultado que mi función.
% function boosted_spinor = BoostBispinor(psi, beta)
% % boostDiracSpinor Boosts a 4-component Dirac spinor to a moving reference frame.
% %
% % Inputs:
% %   psi  - 4x1 complex vector representing the Dirac spinor
% %   beta - 3x1 or 1x3 real vector representing velocity v/c (|beta| < 1)
% %
% % Output:
% %   boosted_spinor - 4x1 complex vector representing the boosted spinor
% 
%     % Reshape inputs
%     psi = psi(:);
%     beta = beta(:);
%     
%     if length(psi) ~= 4
%         error('Input spinor psi must be a 4-element vector.');
%     end
%     if length(beta) ~= 3
%         error('Velocity vector beta must have 3 components.');
%     end
% 
%     beta_sq = sum(beta.^2);
%     beta_norm = sqrt(beta_sq);
% 
%     % Zero velocity edge case
%     if beta_norm == 0
%         boosted_spinor = psi;
%         return;
%     end
% 
%     % Physical boundary check
%     if beta_norm >= 1
%         error('Speed |beta| must be strictly less than 1 (c=1).');
%     end
% 
%     % Calculate Lorentz factor
%     gamma = 1 / sqrt(1 - beta_sq);
% 
%     % Unit velocity vector
%     b_hat = beta / beta_norm;
% 
%     % Pauli matrices
%     sigma_x = [0 1; 1 0];
%     sigma_y = [0 -i; i 0];
%     sigma_z = [1 0; 0 -1];
% 
%     % Projected Pauli term: sigma . b_hat
%     sigma_dot_b = b_hat(1)*sigma_x + b_hat(2)*sigma_y + b_hat(3)*sigma_z;
% 
%     % Build 4x4 matrix: [0, sigma.b; sigma.b, 0]
%     K_dot_b = [zeros(2), sigma_dot_b; ...
%                sigma_dot_b, zeros(2)];
% 
%     % Boost operator S(Lambda) = cosh(chi/2)*I - sinh(chi/2)*K_dot_b
%     cosh_half = sqrt((gamma + 1) / 2);
%     sinh_half = sqrt((gamma - 1) / 2);
% 
%     S = cosh_half * eye(4) - sinh_half * K_dot_b;
% 
%     % Apply boost transform
%     boosted_spinor = S * psi;
% end