% 
% clc; clear; close all;
% 
% % Example usage
% N = 99;
% n = 0:N; 
% % phi = mod(2 * pi * sin(0.1 * n), 2 * pi); % Simulated wrapped phase
% f0 = 1/N;
% sf = 6;
% psi_u = sf*sin(2*pi*f0*n);
% x = exp(1j*psi_u);
% phi = angle(x);
% % unwrapped_phase = unwrap_L2(phi);
% 
% lambda = 0.1; % Regularization parameter
% iterations = 100; % Iteration count
% unwrapped_phase_L1 = unwrap_L1(phi, lambda, iterations);
% 
% figure;
% plot(n, phi, 'r--', 'LineWidth', 1.5); hold on;
% plot(n, unwrapped_phase_L1, 'b', 'LineWidth', 1.5);
% legend('Wrapped Phase', 'Unwrapped Phase (L1)');
% xlabel('Sample Index'); ylabel('Phase (rad)');
% title('1D Phase Unwrapping using L1-norm');
% grid on;
% 
% 
% function psi = unwrap_L1(phi, lambda, iter)
%     N = length(phi);
%     dphi = diff(phi);
%     dpsi = round((dphi + pi) / (2 * pi)) * 2 * pi - pi;
%     g = dphi - dpsi;
% 
%     % Iterative TV minimization (Split Bregman-like)
%     psi = zeros(1, N);
%     for k = 1:iter
%         for i = 2:N
%             psi(i) = psi(i-1) + g(i-1) - lambda * sign(psi(i) - psi(i-1));
%         end
%     end
%     psi = psi + phi(1); % Adjust initial condition
% end



clc; clear; close all;

% Example usage
N = 99;
n = 0:N; 
f0 = 1/N;
sf = 6;
psi_u = sf*sin(2*pi*f0*n);
x = exp(1j*psi_u);  % Simulated complex values
phi = angle(x);     % Wrapped phase

lambda = 0.1; % Regularization parameter
iterations = 100; % Iteration count
unwrapped_phase_L1 = unwrap_L1(phi, lambda, iterations);

figure;
plot(n, phi, 'r--', 'LineWidth', 1.5); hold on;
plot(n, unwrapped_phase_L1, 'b', 'LineWidth', 1.5);
legend('Wrapped Phase', 'Unwrapped Phase (L1)');
xlabel('Sample Index'); ylabel('Phase (rad)');
title('1D Phase Unwrapping using L1-norm');
grid on;

% Function to perform L1-norm phase unwrapping using TV minimization
function psi = unwrap_L1(phi, lambda, iter)
    N = length(phi);
    psi = zeros(1, N);
    dphi = diff(phi);  % Phase differences
    g = dphi;  % Start with the phase gradient as the difference

    % Iterative Total Variation minimization (proximal gradient)
    for k = 1:iter
        for i = 2:N
            % Update the unwrapped phase using a simple form of TV minimization
            psi(i) = psi(i-1) + g(i-1) - lambda * sign(psi(i) - psi(i-1));
        end
    end
    psi = psi + phi(1); % Add the initial phase value to the result
end
