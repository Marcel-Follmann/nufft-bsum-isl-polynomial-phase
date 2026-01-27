close all; clearvars; clc; 

addpath("01_Figures\");
addpath("99_Utilities\utils\");

figureHandler = FigureHandler("savePath","01_Figures\");
figures = {};


nSamples = 20;
nPolyOrder = 3;


pecs = PECS("N",nSamples, "Q",nPolyOrder, "logFile","pecsLog.mat");
[x0, a0] = pecs.getRandomSequence();
pecs.mpsl(x0, a0);

[~, islDb] = pecs.getIsl();

figures{end+1} = figure("Name","ISL Plot"); 
semilogx(islDb);
box on; grid on; 
ylabel("ISL (dB)");
xlabel("Number of MM steps");

for iFigure = 1:numel(figures)
    figureHandler.PrepareEmf(figures{iFigure});
end
clear ans