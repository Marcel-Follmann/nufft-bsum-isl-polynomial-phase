close all; clearvars; clc; 

addpath("99_Utilities\utils\");
addpath("ZoomPlot-MATLAB-main\");

figureHandler = FigureHandler("savePath","01_Figures\");

%% Parameters
nSignal = 20;
paramOptimal = pi/nSignal;
paramMax = pi/2; 
oversampling = 100;

paramStep = paramOptimal / oversampling;
paramVec = (-paramMax + paramStep: ...
	paramStep: paramMax - paramStep).'; % (1 x :)

%% Simulation
% Generate signal
n = 0:nSignal-1;		
signal = exp( 1i*(paramVec.*n.^2) ); % (nSignal x 1)
clear n

% Aperiodic autocorrelation over all parameters
autocorr = zeros(numel(paramVec), 2*nSignal-1);
for iParam = 1:numel(paramVec)
	autocorr(iParam,:) = xcorr( signal(iParam,:) );
end
clear iParam

% ISL = sum |r_xx(k)|^2 with k \neq 0
autocorr(:,nSignal) = 0;
isl = sum( abs(autocorr).^2 ,2);% - abs(autocorr(:, nSignal)).^2;


%% Plot
fig = figure("Name","ISL Minimization");hold all; 

% Lower Bound
plot( [paramVec(1), paramVec(end)]/paramOptimal, ...
	10*log10( nSignal*(nSignal-1) )*ones(2,1), ...
	'--', "DisplayName","Lower Bound");

% ISL
plot(paramVec/paramOptimal, ...
	10*log10(isl), ...
	'.', "DisplayName", "ISL");

xlabel("$$a_2 \cdot \frac{N}{\pi}$$ $$(-)$$");
ylabel("$$\Phi_{ISL}(a_2)$$ $$(dB)$$");

legend("Location","northwest")
grid on; grid minor;

figureHandler.PrepareEmf(fig);
figureHandler.GenerateOutputs();
clear ans