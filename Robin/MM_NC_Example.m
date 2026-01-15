x0 = 2;        % Initial guess
tol = 1e-6;    % Tolerance
max_iter = 100;

[x_vals, f_vals] = mm_minimize_nonconvex(x0, tol, max_iter);

% Plot
figure;
subplot(2,1,1);
plot(0:length(x_vals)-1, x_vals, '-o');
xlabel('Iteration'); ylabel('x value');
title('MM: x progression');

subplot(2,1,2);
plot(0:length(f_vals)-1, f_vals, '-o');
xlabel('Iteration'); ylabel('f(x)');
title('MM: Objective function value');




function [x_vals, f_vals] = mm_minimize_nonconvex(x0, tol, max_iter)

    % Objective Function
    f = @(x) log(1 + exp(x.^2));
    % Gradient of f
    grad_f = @(x) (2*x .* exp(x.^2)) ./ (1 + exp(x.^2));
    
    xk = x0;
    x_vals = xk;
    f_vals = f(xk);
    
    for iter = 1:max_iter
        % Gradient at current point
        grad_fk = grad_f(xk);
        
        % Majorization step: quadratic surrogate
        % We use a Lipschitz constant L for gradient of f (can be adaptive, here fixed)
        L = 10;  % Approximation of the Lipschitz constant for grad_f
        surrogate = @(x) f(xk) + grad_fk*(x - xk) + (L/2)*(x - xk).^2;
        
        % Minimization step: minimize surrogate
        % The minimizer of the quadratic surrogate is:
        x_next = xk - grad_fk / L;
        
        % Save results
        x_vals(end+1) = x_next;
        f_vals(end+1) = f(x_next);
        
        % Check convergence
        if abs(x_next - xk) < tol
            break;
        end
        
        xk = x_next;
    end
end
