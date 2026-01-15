x0 = 1.5;    % Initial guess
tol = 1e-6;
max_iter = 100;

[x_vals, f_vals] = mm_minimize_trig(x0, tol, max_iter);

% Plotting
figure;
subplot(2,1,1);
plot(0:length(x_vals)-1, x_vals, '-o');
xlabel('Iteration'); ylabel('x');
title('MM Optimization of Trigonometric Function');

subplot(2,1,2);
plot(0:length(f_vals)-1, f_vals, '-o');
xlabel('Iteration'); ylabel('f(x)');
title('Objective Function Value');



function [x_vals, f_vals] = mm_minimize_trig(x0, tol, max_iter)
    % Objective function
    f = @(x) sin(x) + 0.5 * x.^2;
    % Gradient of f
    grad_f = @(x) cos(x) + x;
    
    L = 2; % Lipschitz constant for grad_f (since second derivative is bounded by 2)
    
    xk = x0;
    x_vals = xk;
    f_vals = f(xk);
    
    for iter = 1:max_iter
        % Compute gradient
        grad_fk = grad_f(xk);
        
        % MM update
        x_next = xk - grad_fk / L;
        
        % Store values
        x_vals(end+1) = x_next;
        f_vals(end+1) = f(x_next);
        
        % Convergence check
        if abs(x_next - xk) < tol
            break;
        end
        
        xk = x_next;
    end
end
