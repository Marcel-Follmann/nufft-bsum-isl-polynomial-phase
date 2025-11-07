close all; clearvars; clc; 

addpath("99_Utilities\utils\");
figureHandler = FigureHandler("savePath","01_Figures\");
fig = {};

%% Parameters
N = 20;
Q = 3; 

c = rand(Q,1);

%% Simulation
simu = MmApproach(Q, N);

%% Visualization
% fig{end+1} = simu.VisualizeSignal(c);
fig{end+1} = simu.VisualizeAcf(c);
simu.PrintIsl(c);

%% Format Figures
for i = 1:numel(fig)
    figureHandler.PreparePdf(fig{i});
end
figureHandler.GenerateOutputs()
clear ans; 

%% Helper Functions
