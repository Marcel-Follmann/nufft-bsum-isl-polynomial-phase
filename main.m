close all; clearvars; clc; 
addpath("03_Functions\");
addpath("99_Utilities\utils\");
figureHandler = FigureHandler("savePath","01_Figures\");

simulation = SIMULATION();

simulation.run();

fig = simulation.Plot_Function();

figureHandler.PrepareEmf(fig); 
% figureHandler.GenerateOutputs();
