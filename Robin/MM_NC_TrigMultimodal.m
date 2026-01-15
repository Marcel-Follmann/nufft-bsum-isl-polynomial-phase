initial_points = [-3, -1, 0, 1, 3];
tol = 1e-6; max_iter = 100;

for i = 1:length(initial_points)
    x0 = initial_points(i);
    [x_vals, f_vals] = mm_minimize_multimodal(x0, tol, max_iter);
    fprintf('Start at x0=%.2f --> Min at x=%.4f, f(x)=%.4f\n', ...
        x0, x_vals(end), f_vals(end));
end




function [x_vals, f_vals] = mm_minimize_multimodal(x0, tol, max_iter)
    % Objective function
    f = @(x) sin(5*x) + 0.5 * x.^2;
    % Gradient
    grad_f = @(x) 5*cos(5*x) + x;
    
    % Lipschitz constant estimate (max |second derivative| of f(x))
    % Second derivative: -25*sin(5x) + 1 -> bounded by 25+1 = 26
    L = 26;
    
    xk = x0;
    x_vals = xk;
    f_vals = f(xk);
    
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
end
