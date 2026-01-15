%% MM for f(a2) with psi_{knn'} = 2 a2 * k * (n - n')
clear; clc; close all

%% Parameters
N = 128;               % You can change N
max_iter = 1e3;       % Max iterations
tol = 1e-9;          % Tolerance for convergence

%% Sequence design
L = 1;
Q = 3;
M = N / L;

a_coeff = rand(Q, L);
a2_0 = a_coeff(3);          % Initial guess
[x, phi] = create_sequence(a_coeff, M, L);

ri = xcorr(x);
atocrln_i = 10*log10(abs(ri).^2);
ISL_init = 10*log10(sum(abs(ri([1:N-1, N+1:2*N-1])).^2)/N^2);
lag = (1-N:N-1)';

%% Helper Functions

% f(a2): real-valued objective
f_func = @(a2) compute_f(a2, N);

% Gradient of f w.r.t a2
grad_f = @(a2) compute_grad_f(a2, N);

% Estimate Lipschitz constant conservatively
% Use sum of |2k(n - n')| over all terms
LC = estimate_L(N);
% L = 1e4;
ISL_itr = [];
a_itr = [];
%% Check alternate a_coeff values
% a2_alt = [a2_0 - pi/2; a2_0 - pi/4; a2_0; a2_0+pi/4; a2_0+pi/2];
a2_alt = a2_0 + [pi*rand(3,1); -pi*rand(3,1)];
for j = 1:length(a2_alt)
    a2_0 = a2_alt(j);
    %% MM Iterations
    a2_vals = a2_0;
    f_vals = f_func(a2_0);
    a2_curr = a2_0;

    for iter = 1:max_iter
        grad = grad_f(a2_curr);

        % MM Update
        a2_next = a2_curr - grad / LC;

        a2_vals(end+1) = a2_next;
        f_vals(end+1) = f_func(a2_next);

        if abs(a2_next - a2_curr) < tol
            break;
        end
        a2_curr = a2_next;
    end
    fprintf('Optimal a2: %.6f\n', a2_curr);
    fprintf('Optimal f(a2): %.6f\n', f_func(a2_curr));
    %% Plot convergence
    figure;
    subplot(2,1,1);
    plot(0:length(a2_vals)-1, a2_vals, 'ro-','LineWidth',1.5);
    xlabel('Iteration'); ylabel('a_2');
    title('Convergence of a_2');
    
    subplot(2,1,2);
    plot(0:length(f_vals)-1, f_vals, 'bo-','LineWidth',1.5);
    xlabel('Iteration'); ylabel('f(a_2)');
    title('Objective Function Convergence');

    a_coeff_f = a_coeff;
    a_coeff_f(3) = a2_curr;

    [xf, phi] = create_sequence(a_coeff_f, M, L);
    rf = xcorr(xf);
    atocrln_f = 10*log10(abs(rf).^2);
    ISL_opt = 10*log10(sum(abs(rf).^2)/N^2);
    a_itr(end+1) = a2_curr;
    ISL_itr(end+1) = ISL_opt;
end

[~,Indx] = min(ISL_itr);
a_coeff_f(3) = a_itr(Indx);
[xf, phi] = create_sequence(a_coeff_f, M, L);
rf = xcorr(xf);
atocrln_f = 10*log10(abs(rf).^2);
ISL_opt = 10*log10(sum(abs(rf([1:N-1, N+1:2*N-1])).^2)/N^2);
figure; hold on; box on;
plot(lag, atocrln_i, '-.r');
plot(lag, atocrln_f, '-b*');

fprintf('Initial ISL: %f \n', ISL_init);
fprintf('Optimal ISL: %f \n', ISL_opt);

%% --- Function Definitions ---

function f_val = compute_f(a2, N)
    total = 0;
    for k = 1:N-1
        for n = 1:N-k
            for n_prime = 1:N-k
                psi = 2 * a2 * k * (n - n_prime);
                total = total + cos(psi); % only real part contributes
            end
        end
    end
    f_val = total;
end

function grad = compute_grad_f(a2, N)
    total = 0;
    for k = 1:N-1
        for n = 1:N-k
            for n_prime = 1:N-k
                delta = n - n_prime;
                psi = 2 * a2 * k * delta;
                total = total - 2 * k * delta * sin(psi);
            end
        end
    end
    grad = total;
end

function LC = estimate_L(N)
    % Conservative Lipschitz constant estimation
    sum_terms = 0;
    for k = 1:N-1
        for n = 1:N-k
            for n_prime = 1:N-k
                delta = n - n_prime;
                sum_terms = sum_terms + abs(2 * k * delta);
            end
        end
    end
    LC = 10*sum_terms;  % conservative scaling
%     LC = sum_terms;
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