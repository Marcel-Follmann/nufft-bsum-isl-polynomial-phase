close all; clearvars; clc; 

addpath("01_Figures\");
addpath("99_Utilities\utils\");

figureHandler = FigureHandler("savePath","01_Figures\");
figures = {};

%%

nMonteCarlo = 3;
vSamples(:,1) = 10:10:50;
vPolyOrder(:,1) = 3:5;

nSamples = numel(vSamples);
nPolyOrder = numel(vPolyOrder);

timerPECS = zeros(nSamples, nPolyOrder, nMonteCarlo);
timerCdNufft = zeros(nSamples, nPolyOrder, nMonteCarlo);

for iSamples = 1:nSamples
    for iPolyOrder = 1:nPolyOrder
        fprintf("Samples: %d, Order: %d -> ", vSamples(iSamples), vPolyOrder(iPolyOrder));

        % Prepare Methods
        pecs = PECS("N",vSamples(iSamples), "Q",vPolyOrder(iPolyOrder));
        cdmm = CDMM("nr_samples",vSamples(iSamples), "nr_params", vPolyOrder(iPolyOrder), "nr_candidates",1);

        for iMonteCarlo = 1:nMonteCarlo
            [x0, a0] = pecs.getRandomSequence();

            % PECS
            tic; 
            pecs.mpsl(x0, a0);
            timerPECS(iSamples, iPolyOrder, iMonteCarlo) = toc;

            % CDMM
            tic; 
            cdmm.simulate(a0);
            timerCdNufft(iSamples, iPolyOrder, iMonteCarlo) = toc;

            fprintf(".");
        end
        fprintf("\n");
    end
end

meanTimerPECS = mean(timerPECS, 3); 
meanTimerCdNufft = mean(timerCdNufft, 3);

figures{end+1} = figure("Name", "Runtime Comparison Samples"); 
plot(vSamples, meanTimerPECS, "DisplayName","PECS"); hold on; 
plot(vSamples, meanTimerCdNufft, "DisplayName","CDMM");
box on; grid on; 
xlabel("Nr of Samples (-)");
ylabel("Runtime (s)")
legend("Location","best");

figures{end+1} = figure("Name", "Runtime Comparison Order"); 
plot(vPolyOrder, meanTimerPECS, "DisplayName","PECS");hold on; 
plot(vPolyOrder, meanTimerCdNufft, "DisplayName","CD MM");
box on; grid on; 
xlabel("Polynomial Order (-)");
ylabel("Runtime (s)")
legend("Location","best");

for iFigure = 1:numel(figures)
    figureHandler.PrepareEmf(figures{iFigure});
end
clear ans