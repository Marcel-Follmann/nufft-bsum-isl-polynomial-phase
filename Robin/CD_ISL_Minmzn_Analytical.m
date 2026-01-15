% clc; clear; close all;
% rng('shuffle');
% addpath('..\common\')
% 
% %% Sequence Set Design parameters
% N = 128;
% L = 1;
% Q = 3;
% max_itrns = 10;
% M = N/L;
% 
% lags = 1-N:N-1;
% 
% % Seed Sequence design
% a_coeff_o_indx = 3;
% a_coeff_i = rand(Q,L);
% [xi] = create_sequence(a_coeff_i, M, L);
% 
% 
% 
% 
% %% Support functions
% function [x] = create_sequence(a_coeff, M, L)
% theta_li = [];
% m = (1:M).';
% for l = 1:L
%     theta_li(:,l) = a_coeff(1,l)+ a_coeff(2,l).*m + a_coeff(3,l).*m.^2;
% end
% theta_x = theta_li(:);
% x = exp(1j*theta_x);
% end


clc; clear; close all

%% Parameters
N = 64;
maxIter = 100;
tol = 1e-6;

% Initialize psi in 2D form
psi = cell(N-1,1);
for k = 1:N-1
    psi{k} = -pi + 2*pi*rand(N-k, 1);
end

ob_h = zeros(maxIter,1);

%% MM Iterations
for iter = 1:maxIter
    psi_old = psi;
    
    % MM update for each (k,n)
    for k = 1:N-1
        for n = 1:N-k
            psi{k}(n) = psi_old{k}(n) + sin(psi_old{k}(n));
        end
    end
    
    % Objective function computation: double sum
    obj = 0;
    for k = 1:N-1
        for n = 1:N-k
            obj = obj + cos(psi{k}(n));
        end
    end
    ob_h(iter) = obj;
    
    % Convergence check
    max_diff = 0;
    for k = 1:N-1
        diff_kn = abs(psi{k} - psi_old{k});
        max_diff = max(max_diff, max(diff_kn));
    end
    if max_diff < tol
        fprintf('Converged at iteration %d with max change %.6f\n', iter, max_diff);
        ob_h = ob_h(1:iter);
        break
    end
end

fprintf('Final objective value: %.6f\n', obj);

%% Plot convergence
figure;
plot(1:length(ob_h), ob_h, '-o', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Objective Value');
title('Convergence of Double Sum Objective');
grid on;

