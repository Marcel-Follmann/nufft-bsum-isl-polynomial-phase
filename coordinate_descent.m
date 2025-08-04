% coordinate_descent_surface_clean.m
% Zwei 3-D-Höhenplots mit Coordinate-Descent-Pfad



clear; clc; close all;

figureHandler = FigureHandler("savePath","01_Figures\");

% Parameter -----------------------------------------------------------
x0      = [-3; 2];                     % Startpunkt
A1 = [5 0; 0 1];
A2 = [5 2; 2 1];


% Plot ----------------------------------------------------------------
fig1 = Visualize(A1, x0, "Coordinate Descent uncoupled");
fig2 = Visualize(A2, x0, "Coordinate Descent coupled");

figureHandler.PrepareEmf(fig1, "width", 15.5);
figureHandler.PrepareEmf(fig2, "width", 15.5);
figureHandler.GenerateOutputs();

function fig = Visualize(A, x0, title)
	fig = figure("Name",title);
	[x, y, z] = GetEllipse(A);
	contour(x, y, z, "DisplayName","$$f(x)$$"); hold on; 
	shading interp
	
	[x, y, z] = CoordinateDescent(A, x0, 10);
	
	
	plot3(x(end),y(end),z(end),'s', "DisplayName","$$x^*$$");
	plot3(x(1),  y(1),  z(1),  's', "DisplayName","$$x_0$$");
	plot3(x, y, z, '-k', "DisplayName","$$x_k$$");
	grid on; grid minor;
	legend("Location","northeast")
	
	xlabel("x coordinate");
	ylabel("y coordinate");
end


function [x, y, z] = GetEllipse(A)
	[x,y] = meshgrid(-4:0.1:4);
	z = 0.5*(A(1,1)*x.^2 + 2*A(1,2)*x.*y + A(2,2)*y.^2);
end

function [x, y, z] = CoordinateDescent(A, x0, iter_max)
	x = x0;
	path = x;
	z_tmp = Inf;
	i = 0;
	tol = 1e-3;

	while (abs(z_tmp) >= tol) && (i < iter_max)
		i = i+1;
		% x-Update
		x(1)  = -(A(1,2)*x(2))/A(1,1); 
		path  = [path x]; %#ok<AGROW>

        % y-Update
        x(2)  = -(A(2,1)*x(1))/A(2,2);
        path  = [path x]; %#ok<AGROW>
		z_tmp = 0.5*sum((A*x).*x,1);
	end
	
	x = path(1,:);
	y = path(2,:);
	z = 0.5*sum((A*path).*path,1);
end