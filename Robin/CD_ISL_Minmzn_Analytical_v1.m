clc; clear; close all

%% Parameters
N = 128;
L = 1;
Q = 3;
M = N/L;

maxIter = 100;
tol = 1e-6;

a_coeff = rand(Q,L);
[x, phi] = create_sequence(a_coeff, M, L);

%% MM Iterations
ob_h = zeros(maxIter,1);
psi_old = 0;

for iter = 1:maxIter  
    % MM update for each (k,n) & Objective function computation: double sum
    obj = 0;
    for k = 1:N-1
        for n = 1:N-k
            for np = 1:N-k
                if np ~= n
                    psi_knnp = 2*a_coeff(3,L)*k*(n-np);
                    obj = obj + cos(psi_knnp);
                else
                    continue
                end
            end
        end
    end
    ob_h(iter) = obj;
    
    % Convergence check
    max_diff = 0;
    for k = 1:N-1
        diff_kn = abs(psi_knnp - psi_old);
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

%% Support functions
function [x, phi_x] = create_sequence(a_coeff, M, L)
phi_i = [];
m = (1:M).';
for l = 1:L
    phi_i(:,l) = a_coeff(1,l)+ a_coeff(2,l).*m + a_coeff(3,l).*m.^2;
end
phi_x = phi_i(:);
x = exp(1j*phi_x);
end