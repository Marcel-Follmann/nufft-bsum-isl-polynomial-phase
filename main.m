close all; clearvars; clc; 
addpath("99_Utilities\utils\");

figureHandler = FigureHandler("savePath","01_Figures\");
figures = {};

%% Parameters
maxIterations = 20;
N = 20;
Q = 4;
stepImprovementGoal_dB = 0.1;

a0 = randn(Q,1);

%% Coordinate Descent Simulation
simulation = Simulation( ...
	"N",N, ...
	"Q",Q ...
	);

a(:,1) = a0;
isl(1) = simulation.ComputeIsl_dB(a(:,1));
for i = 1:maxIterations

	a(:,i+1) = simulation.CoordinateDescent(a(:,i), ...
		"fftOversampling",10 ...
		);

	isl(i+1) = simulation.ComputeIsl_dB(a(:,i));
end

%% Plot
figures{end+1} = figure("Name", "Convergence");
plot(0:maxIterations,10*log10(isl), ...
	'-', "DisplayName", "ISL");
xlabel("CD Iterations ($-$)");
ylabel("ISL ($dB$)");
grid on;

fprintf("Initial ISL: %.2e\n", isl(1));
fprintf("CD ISL: %.2e\n", isl(end));

%% Format Figures
for iFigure = 1:numel(figures)
	figureHandler.PreparePdf(figures{iFigure}, ...
		"width",		32,			...
		"height",		16,			...
		"font",			"Arial",	...
		"fontSize",		20,			...
		"lineWidth",	3,			...
		"markerSize",	8			...
		);
end
figureHandler.GenerateOutputs();
clear ans

%% Helper Functions

