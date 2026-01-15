close all; clearvars; clc; 

addpath("99_Utilities\utils\");
figureHandler = FigureHandler("savePath","01_Figures\");
fig = {};

%% Parameters
N = 20;
Q = 4; 

nStartpoints = 15;
itersCD = 20; 
itersMM = 100;

cStart = rand(Q,nStartpoints);



%% Simulation
simu = MmApproach(Q, N);

isl = NaN * ones(itersCD,itersMM,nStartpoints);
for iStartpoint = 1:nStartpoints
    c = cStart(:,iStartpoint);
    isl(1,1,iStartpoint) = 10*log10(simu.DeriveIsl(c));

    for iCD = 1:itersCD
        for p = 2:Q-1
            for iMM = 1:itersMM
                [c, go2NextCoordinate] = simu.MmStep(p, c, "Tol",1e-6);
                isl(iCD, iMM, iStartpoint) = 10*log10(simu.DeriveIsl(c)); 

                if iMM > 2
                    if abs(isl(iCD, iMM-1:iMM, iStartpoint)) <= 1e-1
                        break;
                    end
                end
            end
        end
        islEnd(iStartpoint) = isl(iCD, iMM, iStartpoint);
    end
    cStart(:,iStartpoint) = c;
end

%% Visualization
% Logging
for iStartpoint = 1:nStartpoints
    fprintf("Run %2d, ISL: %.2f dB -> %.2f dB (%5.2f dB)\n", ...
        iStartpoint, ...
        isl(1,1, iStartpoint), ...
        islEnd(iStartpoint), ...
        isl(1,1, iStartpoint) - islEnd(iStartpoint));
end
islWorst = max(isl(1,1,:));
[islBest,iBest] = min(islEnd);
fprintf("Isl worst %.2f dB, best %.2f dB ( %.2f dB)\n", islWorst, islBest, islWorst - islBest);

% Plot ISL
% fig{end+1} = figure("Name","ISL History"); hold on
for iStartpoint = 1:nStartpoints
    linetype = '-.';
    if iStartpoint == iBest
        linetype = '-';
    end
    [X,Y] = meshgrid(1:itersCD, 1:itersMM);
    plot3(X, Y, isl(:,:,iStartpoint));%, linetype, ...
        % "DisplayName","Run " + num2str(iStartpoint));
end
grid on; grid minor; 
ylabel("ISL (dB)");
xlabel("Steps (-)");
legend("Location","northeastoutside");
% fig{end+1} = simu.VisualizeSignal(c);
% fig{end+1} = simu.VisualizeAcf(c);
% simu.PrintIsl(c);

%% Format Figures
for i = 1:numel(fig)
    figureHandler.PrepareEmf(fig{i});
end
figureHandler.GenerateOutputs()
clear ans; 

%% Helper Functions
