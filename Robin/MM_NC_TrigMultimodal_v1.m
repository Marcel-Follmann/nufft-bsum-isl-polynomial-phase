clc; clear; close all;

% Objective and gradient
flg_obj = 4;
switch flg_obj
    case 1
        % First Objective
        f = @(x) sin(5*x) + 0.5 * x.^2;
        grad_f = @(x) 5*cos(5*x) + x;
        % MM parameter
        L = 26;              % Lipschitz constant estimate
    case 2
        % Second Objective
        f = @(x) sin(5*x);
        grad_f = @(x) 5*cos(5*x);
        % MM parameter
        L = 26;              % Lipschitz constant estimate
    case 3
        % Third Objective
        f = @(x) cos(x);
        grad_f = @(x) -sin(x);
        % MM parameter
        L = 26;              % Lipschitz constant estimate
    case 4
        % Third Objective
        f = @(x) cos(x) + cos(4*x) + cos(15*x);
        grad_f = @(x) -sin(x) - 4*sin(4*x) - 15*sin(15*x);
        % MM parameter
%         L = +242;              % Lipschitz constant estimate
        L = +150;              % erroneous Lipschitz constant estimate
end


x0 = 10*rand;            % Initial point (change to test different starts)
tol = 1e-6;
max_iter = 1000;

% Initialize
xk = x0;
x_vals = xk;
f_vals = f(xk);

% MM iterations
for iter = 1:max_iter
    grad_fk = grad_f(xk);
    x_next = xk - grad_fk / L;
    
    x_vals(end+1) = x_next;
    f_vals(end+1) = f(x_next);
    
    if abs(x_next - xk) < tol
        break;
    end
    xk = x_next;
end

%% Plot the Objective Function and Optimization Path
x_plot = linspace(-5*x0, 5*x0, 10000);
f_plot = f(x_plot);

figure;

% Plot Objective function
subplot(2,1,1);
plot(x_plot, f_plot, 'b-', 'LineWidth', 1.5); hold on;
plot(x_vals, f(x_vals), 'ro-', 'LineWidth', 1.5);
plot(x_vals(end), f(x_vals(end)), 'ks', 'MarkerSize',10,'MarkerFaceColor','g');
xlabel('x'); ylabel('f(x)');
title('Objective Function with MM Optimization Path');
legend('f(x)', 'MM Path', 'Optimal Point','Location','best');
grid on;

% Plot convergence of objective value
subplot(2,1,2);
plot(0:length(f_vals)-1, f_vals, 'ko-', 'LineWidth', 1.5);
xlabel('Iteration'); ylabel('f(x)');
title('Convergence of Objective Function Value');
grid on;

% Display optimal point
fprintf('Optimal x found: %.4f\n', x_vals(end));
fprintf('Optimal f(x): %.4f\n', f_vals(end));
