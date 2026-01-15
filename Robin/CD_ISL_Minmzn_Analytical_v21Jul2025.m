clc; clear; close all;

%% Parameters
N = 16;
L = 1;
Q = 3;
M = N / L;

maxIter = 10;
tol = 1e-6;

a_coeff = rand(Q, L);
[x, phi] = create_sequence(a_coeff, M, L);

ri = xcorr(x);
atocrln_i = 10*log10(abs(ri).^2);
ISL_init = 10*log10(sum(abs(ri).^2)/N^2);

%% MM iterations
a2 = a_coeff(3);  % Extract initial a2
obj_itr = zeros(maxIter, 1);
g = @(x) cos(x);
for iter = 1:maxIter
    psi_knnp_n = 0;
    obj = 0;
    % MM Update step for each (k, n, n')
    for k = 1:N-1
        for n = 1:N-k
            for np = 1:N-k
                if n < np
                    psi_knnp = 2 * a2 * k * (n - np);
%                     psi_knnp_n = psi_knnp_n + psi_knnp + sin(psi_knnp);
                    obj = obj + 2*g(psi_knnp);
                else
                    continue
                end
            end
        end
    end
    obj_itr(iter) = obj;
end

%% Plot normalized ISL convergence
figure;
plot(obj_itr, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Iteration');
ylabel('Obj');
title('Convergence of Obj via MM');

%% Support functions
function [x, phi_x] = create_sequence(a_coeff, M, L)
phi_i = [];
m = (1:M).';
for l = 1:L
    phi_i(:,l) = a_coeff(1,l) + a_coeff(2,l).*m + a_coeff(3,l).*m.^2;
end
phi_x = phi_i(:);
x = exp(1j*phi_x);
end


%     if iter > 1 && abs(normalized_ISL(iter) - normalized_ISL(iter-1)) < tol
%     if iter > 1 && norm(abs(a2_i) - abs(a2)) < tol
%         fprintf('Converged at iteration %d\n', iter);
%         normalized_ISL = normalized_ISL(1:iter);
%         break;
%     end