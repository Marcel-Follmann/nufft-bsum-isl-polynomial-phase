close all; clearvars; clc; 

addpath("99_Utilities\utils\");
addpath("01_Figures\");
figureHandler = FigureHandler("savePath","01_Figures\");

%%

nSamples = 2^8;
nNeighbours = 6;
omega = 5;
fZeroPadding = 10;

X = zeros(nSamples,1);
omegaVec = 0:nSamples-1;

idx = find(abs(omegaVec) == omega);
X(idx) = 1/nSamples;

dOmega = omegaVec(2) - omegaVec(1);
T = 1/dOmega;
dt = 1/nSamples;
t = 0:dt:T-dt; 

x = fftshift(ifft(X))*nSamples;
xZeroPadding = fftshift(ifft(X,fZeroPadding*nSamples))*nSamples*fZeroPadding;
tZeroPadding = 0:dt/fZeroPadding:T-dt/fZeroPadding; 

tw = tZeroPadding;%T*rand(nSamples,1);
nufft = NUFFT("K",2*nSamples, "x",X);
nufft.prepareInterpolator("J",6);
Xw = nufft.interpolate(tw*2*pi);
for i = 1:numel(tw)
    
    e(i) = nufft.error(tw(i)*2*pi);
end

% Xw = interp_fft_gaussian(x, tw, nNeighbours, 1);

nexttile;
plot(omegaVec, X, 'o');
ylim([0,1]);
nexttile; 
plot(t, real(x), '.'); hold on;
plot(tZeroPadding, real(xZeroPadding), '--');
plot(tw, real(Xw), 'o');
nexttile; hold on; 
plot(tZeroPadding, abs(real(Xw.') - real(xZeroPadding)));
plot(tZeroPadding, e);
yLim(:,1) = ylim;
% plot([t;t], yLim*ones(size(t)), 'k')

