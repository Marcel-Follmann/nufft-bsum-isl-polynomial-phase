close all; clearvars; clc; 

addpath("99_Utilities\utils\");
addpath("01_Figures\");
figureHandler = FigureHandler("savePath","01_Figures\");

%%

nSamples = 2^9;
nNeighbours = 6;
omega = 5;
fZeroPadding = 4;

X = zeros(nSamples,1);
omegaVec = -nSamples/2:nSamples/2-1;

idx = find(abs(omegaVec) == omega);
X(idx) = 1/sqrt(8);

idx = find(abs(omegaVec) == 5*omega);
X(idx) = 1/sqrt(8);

idx = find(abs(omegaVec) == 15*omega);
X(idx) = 1/sqrt(8);
idx = find(abs(omegaVec) == 25*omega);
X(idx) = 1/sqrt(8);

dOmega = omegaVec(2) - omegaVec(1);
T = 1/dOmega;
dt = 1/nSamples;
t = 0:dt:T-dt; 

x = real(fft(ifftshift(X)))/sqrt(nSamples);

pad = (fZeroPadding-1)*nSamples/2;
Xpadded = [zeros(pad, 1); X; zeros(pad, 1)];
xZeroPadding = real(fft(ifftshift(Xpadded)))/sqrt(nSamples);
tZeroPadding = 0:dt/fZeroPadding:T-dt/fZeroPadding; 

tw = tZeroPadding;%T*rand(nSamples,1);
nufft = NUFFT("K",fZeroPadding*nSamples, "x",X);
nufft.prepareInterpolator("J",nNeighbours);
Xw = real(nufft.interpolate(tw*2*pi));

e = (nufft.error(tw*2*pi));


% Xw = interp_fft_gaussian(x, tw, nNeighbours, 1);

nexttile;
plot(omegaVec, X, 'o', "DisplayName","Spectral Function");
ylim([0,1]);
legend("Location","best")
nexttile; 
plot(t, x, '.', "DisplayName","FFT"); hold on;
plot(tZeroPadding, xZeroPadding, '--', "DisplayName","Zero Padded Version");
plot(tw, real(Xw), 'o', "DisplayName","Interpolated Version");
legend("Location","best")
nexttile; hold on; 
plot(tZeroPadding, 10*log10(abs(real(Xw) - real(xZeroPadding))), "DisplayName","Actual Error");
plot(tZeroPadding, 10*log10(e), "DisplayName","Empirical Worst Case");
plot([tZeroPadding(1), tZeroPadding(end)], 10*log10(10*eps*max(x)*[1,1]), "DisplayName","Machine Precision");
legend("Location","best")
% yLim(:,1) = ylim;
% plot([t;t], yLim*ones(size(t)), 'k')

