close all; clearvars; clc; 
addpath("99_Utilities\utils\");

figureHandler = FigureHandler("savePath","01_Figures\");
figures = {};

%% Parameters
maxIterations = 10;
N = 20;
Q = 4;

a0 = [-0.920274947429725;1.179610379138849;1.508569156483106;1.062646144475084];%randn(Q,1);

%% Coordinate Descent Simulation
simulation = Simulation( ...
	"N",N, ...
	"Q",Q ...
	);

a(:,1) = a0;
isl_dB(1) = simulation.ComputeIsl_dB(a(:,1));
for i = 1:maxIterations

	a(:,i+1) = simulation.CoordinateDescent(a(:,i), ...
		"fftOversampling",1 ...
		);

	isl_dB(i+1) = simulation.ComputeIsl_dB(a(:,i));
end

%% Plot
figures{end+1} = figure("Name", "Convergence");
semilogx(0:maxIterations,isl_dB, ...
	'-o', "DisplayName", "ISL");
xlabel("CD Iterations ($-$)");
ylabel("ISL ($dB$)");
grid on;

fprintf("Initial ISL: %.2e\n", isl_dB(1));
fprintf("CD ISL: %.2e\n", isl_dB(end));

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

