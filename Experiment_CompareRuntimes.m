close all; clearvars; clc; 

addpath("01_Figures\");
addpath("99_Utilities\utils\");

figureHandler = FigureHandler("savePath","01_Figures\");
figures = {};

%%

nMonteCarlo = 100;
vSamples(:,1) = 20:20:100;
vPolyOrder(:,1) = 3:5;

nSamples = numel(vSamples);
nPolyOrder = numel(vPolyOrder);

timerPECS = zeros(nSamples, nPolyOrder, nMonteCarlo);
for iSamples = 1:nSamples
    for iPolyOrder = 1:nPolyOrder
        pecs = PECS("N",iSamples, "Q",iPolyOrder);
        for iMonteCarlo = 1:nMonteCarlo
            [x0, a0] = pecs.getRandomSequence();
            tic; 
            pecs.mpsl(x0, a0);
            timerPECS(iSamples, iPolyOrder, iMonteCarlo) = toc;
        end
    end
end

meanTimerPECS = mean(timerPECS, 3); 

figures{end+1} = figure("Name", "Runtime Comparison");

plot(vPolyOrder, meanTimerPECS);
xlabel("Polynomial Order (-)");
ylabel("Runtime (s)")
legend("Nr. of Samples: " + string(vSamples), "Location","best");

for iFigure = 1:numel(figures)
    figureHandler.PrepareEmf(figures{iFigure});
end
clear ans