close all; clearvars; clc; 

addpath("01_Figures\");
addpath("99_Utilities\utils\");

figureHandler = FigureHandler("savePath","01_Figures\");
figures = {};

%%

nMonteCarlo = 1;
vSamples(:,1) = 10:5:20;
vPolyOrder(:,1) = 3:5;

nSamples = numel(vSamples);
nPolyOrder = numel(vPolyOrder);

timerPECS = zeros(nSamples, nPolyOrder, nMonteCarlo);
for iSamples = 1:nSamples
    for iPolyOrder = 1:nPolyOrder
        fprintf("Samples: %d, Order: %d -> ", vSamples(iSamples), vPolyOrder(iPolyOrder));
        pecs = PECS("N",vSamples(iSamples), "Q",vPolyOrder(iPolyOrder));
        for iMonteCarlo = 1:nMonteCarlo
            [x0, a0] = pecs.getRandomSequence();
            tic; 
            pecs.mpsl(x0, a0);
            timerPECS(iSamples, iPolyOrder, iMonteCarlo) = toc;
            fprintf(".");
        end
        fprintf("\n");
    end
end

meanTimerPECS = mean(timerPECS, 3); 


figures{end+1} = figure("Name", "Runtime Comparison Samples");
plot(vSamples, meanTimerPECS);
xlabel("Nr of Samples (-)");
ylabel("Runtime (s)")
legend("Order: " + string(vPolyOrder), "Location","best");

figures{end+1} = figure("Name", "Runtime Comparison Order");
plot(vPolyOrder, meanTimerPECS);
xlabel("Polynomial Order (-)");
ylabel("Runtime (s)")
legend("Nr of Samples: " + string(vSamples), "Location","best");

for iFigure = 1:numel(figures)
    figureHandler.PrepareEmf(figures{iFigure});
end
clear ans