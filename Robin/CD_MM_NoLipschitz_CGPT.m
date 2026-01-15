%% MM with Cosine Majorizer (no Lipschitz constant)
clear; clc; close all

%% Parameters
N = 128;
max_iter = 1e3;
tol = 1e-6;

%% Sequence design
L = 1;
Q = 3;
M = N / L;

a_coeff = rand(Q, L);
a2_0 = a_coeff(3);          % Initial guess
[x, phi] = create_sequence(a_coeff, M, L);

ri = xcorr(x);
atocrln_i = 10*log10(abs(ri).^2./N^2);
ISL_init = 10*log10(sum(abs(ri).^2./N^2));
lag = 1-N:N-1;

%% Initialize
a2_vals = a2_0;
f_vals = compute_f(a2_0, N);
a2_curr = a2_0;

for iter = 1:max_iter
    B = 0;
    C = 0;
    
    for k = 1:N-1
        for n = 1:N-k
            for n_prime = 1:N-k
                d = 2 * k * (n - n_prime);
                psi_curr = d * a2_curr;
                B = B - d * sin(psi_curr);
                C = C + d^2;
            end
        end
    end
    
    % MM update (minimize the surrogate)
    a2_next = a2_curr - B / C;
    
    % Store
    a2_vals(end+1) = a2_next;
    f_vals(end+1) = compute_f(a2_next, N);
    
    if abs(a2_next - a2_curr) < tol
        break;
    end
    
    a2_curr = a2_next;
end

fprintf('Optimal a2: %.6f\n', a2_curr);
fprintf('Optimal f(a2): %.6f\n', f_vals(end));

%% Plot convergence
figure;
subplot(2,1,1);
plot(0:length(a2_vals)-1, a2_vals, 'ro-','LineWidth',1.5);
xlabel('Iteration'); ylabel('a_2');
title('Convergence of a_2 (Majorizer without Lipschitz)');

subplot(2,1,2);
plot(0:length(f_vals)-1, f_vals, 'bo-','LineWidth',1.5);
xlabel('Iteration'); ylabel('f(a_2)');
title('Objective Function Convergence');


a_coeff_f = a_coeff;
a_coeff_f(3) = a2_curr;
[xf, phi] = create_sequence(a_coeff_f, M, L);
rf = xcorr(xf);
atocrln_f = 10*log10(abs(rf).^2./N^2);
ISL_opt = 10*log10(sum(abs(rf).^2./N^2));

figure; hold on; box on;
plot(lag, atocrln_i, '-.r');
plot(lag, atocrln_f, '-b*');

fprintf('Initial ISL: %f \n', ISL_init);
fprintf('Optimal ISL: %f \n', ISL_opt);
%% Function to compute f(a2)
function f_val = compute_f(a2, N)
    total = 0;
    for k = 1:N-1
        for n = 1:N-k
            for n_prime = 1:N-k
                psi = 2 * a2 * k * (n - n_prime);
                total = total + cos(psi);
            end
        end
    end
    f_val = total;
end


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