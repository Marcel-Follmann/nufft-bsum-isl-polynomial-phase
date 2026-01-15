%% MM Algorithm with Surrogate Plotting and Monotonicity Check
clear; clc; close all

%% Parameters
N = 5;
max_iter = 50;
tol = 1e-6;
a2_0 = 0.5;

%% Initialization
a2_vals = a2_0;
f_vals = compute_f(a2_0, N);
a2_curr = a2_0;
non_monotonic_count = 0;

%% Plot range
a2_range = linspace(-2, 2, 400);

%% Main MM Loop
for iter = 1:max_iter
    B = 0;
    C = 0;
    d_list = [];

    for k = 1:N-1
        for n = 1:N-k
            for n_prime = 1:N-k
                d = 2 * k * (n - n_prime);
                d_list(end+1) = d;
                psi_curr = d * a2_curr;
                B = B - d * sin(psi_curr);
                C = C + d^2;
            end
        end
    end

    % MM update
    a2_next = a2_curr - B / C;

    % Evaluate f and check monotonicity
    f_next = compute_f(a2_next, N);
    if f_next > f_vals(end)
        fprintf('⚠️ Warning: Non-monotonic increase detected at iteration %d\n', iter);
        non_monotonic_count = non_monotonic_count + 1;
    end

    % Store values
    a2_vals(end+1) = a2_next;
    f_vals(end+1) = f_next;

    % Plot original function and surrogate
    f_plot = arrayfun(@(a2) compute_f(a2, N), a2_range);
    g_plot = zeros(size(a2_range));
    for i = 1:length(a2_range)
        delta = a2_range(i) - a2_curr;
        surrogate = 0;
        for d = d_list
            psi_curr = d * a2_curr;
            surrogate = surrogate + ...
                cos(psi_curr) - d * sin(psi_curr) * delta - 0.5 * d^2 * delta^2;
        end
        g_plot(i) = surrogate;
    end

    figure(1); clf;
    plot(a2_range, f_plot, 'b', 'LineWidth', 2); hold on;
    plot(a2_range, g_plot, 'r--', 'LineWidth', 2);
    plot(a2_curr, compute_f(a2_curr, N), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    plot(a2_next, f_next, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    legend('Original f(a_2)', 'Majorizer g(a_2)', 'Current point', 'Next point','Location','southeast');
    title(sprintf('Iteration %d', iter));
    xlabel('a_2'); ylabel('Function value');
    drawnow;

    if abs(a2_next - a2_curr) < tol
        break;
    end
    a2_curr = a2_next;
end

fprintf('\n✅ Final optimal a2 = %.6f\n', a2_curr);
fprintf('✅ Final objective f(a2) = %.6f\n', f_vals(end));
fprintf('🔁 Total non-monotonic occurrences: %d\n', non_monotonic_count);

%% Final convergence plots
figure;
subplot(2,1,1);
plot(0:length(a2_vals)-1, a2_vals, 'ro-','LineWidth',1.5);
xlabel('Iteration'); ylabel('a_2');
title('Convergence of a_2');

subplot(2,1,2);
plot(0:length(f_vals)-1, f_vals, 'bo-','LineWidth',1.5);
xlabel('Iteration'); ylabel('f(a_2)');
title('Objective Function f(a_2) per Iteration');

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
