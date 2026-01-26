close all; clearvars; clc; 
addpath("99_Utilities\utils\");

figureHandler = FigureHandler("savePath","01_Figures\");
figures = {};

%% Parameters
maxIterations = 1;
N = 100;
Q = 4;

a0(:,1) = randn(Q,1);
% a0 = [randn(3,1);zeros(Q-3,1)  ]; % Tends to produce better ISLs

%% Coordinate Descent Simulation
tic;
simulation = Simulation( ...
	"N",N, ...
	"Q",Q ...
	);

a(:,1) = a0;
isl_dB(1) = simulation.ComputeIsl_dB(a(:,1));
for i = 1:maxIterations
	a(:,i+1) = simulation.CoordinateDescent(a(:,i));
	isl_dB(i+1) = simulation.ComputeIsl_dB(a(:,i+1));
end
toc;
lags = -(N-1):N-1;
y0 = xcorr(Signal(a0, N))/N;
y = xcorr(Signal(a(:,end),N))/N;

%% Plot
figures{end+1} = figure("Name", "Convergence");
plot(0:maxIterations,isl_dB, ...
	'-o', "DisplayName", "ISL");
xlabel("CD Iterations ($-$)");
ylabel("ISL ($dB$)");
grid on;

figures{end+1} = figure("Name", "Autocorrelation"); hold on; 
plot(lags, 20*log10(abs(y0)),'-', "DisplayName","Initial");
plot(lags, 20*log10(abs(y)), "DisplayName","Optimized");


fprintf("Initial ISL: %.2f dB\n", isl_dB(1));
fprintf("CD ISL: %.2f dB\n", isl_dB(end));
fprintf("Improvement: %.2f dB\n", isl_dB(1) - isl_dB(end));
legend("Location","best")

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
function s = Signal(a, N)
    n(:,1) = 0:N-1; 
    phase = zeros(N,1);
    for q = 1:numel(a)
        phase = phase + a(q).*n.^(q-1);
    end
    s = exp(1i*phase);
end