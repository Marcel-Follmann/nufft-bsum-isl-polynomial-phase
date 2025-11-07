close all; clearvars; clc; 

addpath("99_Utilities\utils\");
figureHandler = FigureHandler("savePath","01_Figures\");
fig = {};

%% Parameters
N = 20;
Q = 4; 

c = rand(Q,1);

itersCD = 5; 
itersMM = 10;

%% Simulation
simu = MmApproach(Q, N);

islInit = simu.SimulateIsl(c);
i = 1;
for iCD = 1:itersCD
    for p = 2:Q-1
        for iMM = 1:itersMM
            c(:,iMM, iCD) = simu.MmStep(p, c);
            islHist(i) = simu.DeriveIsl(c(:,iMM, iCD));
            i = i+1;
        end
    end
end
islEnd = simu.SimulateIsl(c(:,end, end));

fprintf("ISL Start:   %.2f dB\n", 10*log10(islInit));
fprintf("ISL End:     %.2f dB\n", 10*log10(islEnd));
fprintf("Improvement: %.2f dB\n", 10*log10(islInit) - 10*log10(islEnd));

plot(10*log10(islHist));


%% Visualization
% fig{end+1} = simu.VisualizeSignal(c);
% fig{end+1} = simu.VisualizeAcf(c);
% simu.PrintIsl(c);

%% Format Figures
for i = 1:numel(fig)
    figureHandler.PreparePdf(fig{i});
end
figureHandler.GenerateOutputs()
clear ans; 

%% Helper Functions
