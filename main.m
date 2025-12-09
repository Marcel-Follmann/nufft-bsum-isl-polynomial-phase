close all; clearvars; clc; 
addpath("99_Utilities\utils\");

figureHandler = FigureHandler("savePath","01_Figures\");
figures = {};

%% Parameters
maxIterations = 3;
N = 4;
Q = 4;

a0 = randn(Q,1);
% a0 = [randn(3,1);zeros(Q-3,1)  ]; % Tends to produce better ISLs

%% Coordinate Descent Simulation
simulation = Simulation( ...
	"N",N, ...
	"Q",Q ...
	);

a(:,1) = a0;
isl_dB(1) = simulation.ComputeIsl_dB(a(:,1));
for i = 1:maxIterations

	a(:,i+1) = simulation.CoordinateDescent(a(:,i));

	isl_dB(i+1) = simulation.ComputeIsl_dB(a(:,i));
end

%% Plot
figures{end+1} = figure("Name", "Convergence");
semilogx(0:maxIterations,isl_dB, ...
	'-o', "DisplayName", "ISL");
xlabel("CD Iterations ($-$)");
ylabel("ISL ($dB$)");
grid on;

fprintf("Initial ISL: %.2f dB\n", isl_dB(1));
fprintf("CD ISL: %.2f dB\n", isl_dB(end));
fprintf("Improvement: %.2f dB\n", isl_dB(1) - isl_dB(end));

%% Format Figures
for iFigure = 1:numel(figures)
	figureHandler.PreparePdf(figures{iFigure}, ...
		"width",		32,			...
		"height",		16,			...
		"font",			"Arial",	...
		"fontSize",		20,			...
		"lineWidth",	3,			...
		"markerSize",	4			...
		);
end
figureHandler.GenerateOutputs();
clear ans

%% Helper Functions

