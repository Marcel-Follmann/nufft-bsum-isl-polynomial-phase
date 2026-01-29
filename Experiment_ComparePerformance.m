close all; clearvars; clc; 

addpath("01_Figures\");
addpath("02_Data\");
addpath("99_Utilities\utils\");

figureHandler = FigureHandler("savePath","01_Figures\");
figures = {};


nSamples = 100;
nPolyOrder = 4;


pecs = PECS("N",nSamples, "Q",nPolyOrder, "logFile","pecsLog.mat");
cdmm = CDMM("nr_samples",nSamples, "nr_params",nPolyOrder, "save_data",true, "file_name","cdmmLog", "max_mm_iters",1e3, "nr_candidates",5);

[x0, a0] = pecs.getRandomSequence();
a0 = mod(a0, 2*pi);
pecs.mpsl(x0, a0);
cdmm.simulate(a0);

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