close all; clearvars; clc; 

addpath("01_Figures\");
addpath("02_Data\");
addpath("99_Utilities\utils\");

figureHandler = FigureHandler("savePath","01_Figures\");
figures = {};

cdmm = matfile("02_Data\cdmmLog.mat");
params = cdmm.params;
nufft = cdmm.nufft;

% Reconstruct a over the Optimization
a = cdmm.a(:,1) .* ones(size(cdmm.a,1), size(cdmm.a_mm,2), size(cdmm.a_mm,1));
a(3,:,:) = permute(cdmm.a_mm(:,:,1), [3, 2, 1]);
a_mm = permute(a(3,:,:), [2,3,1]);

isl_dB = get_isl(a, params);

figures{end+1} = figure("Name", "ISL in Coordinate 1"); hold on; 
plot(nufft{1,1}.f*2, 10*log10(nufft{1,1}.X), '.');
plot(a_mm/pi, isl_dB);
box on; grid on; 

figures{end+1} = figure("Name", "Best ISL in Coordinate 1"); 
[~,k] = min(isl_dB(end,:));
semilogx(isl_dB);
box on; grid on; 

for iFigure = 1:numel(figures)
    figureHandler.PrepareEmf(figures{iFigure});
end
clear ans

function isl_dB = get_isl(a, params)
    samples(:,1) = 1:params.nr_samples;
    V = samples .^ (0:params.nr_params-1); % Vandermonde
    phase = pagemtimes(V,a);
    
    x = exp(1i*phase);
    
    F = fft(x, 2*params.nr_samples, 1); % 2N for Aperiodic ACF
    rCirc = ifft(abs(F).^2, [],1);
    rPos = rCirc(2:params.nr_samples,:,:); % lags 1, .., nSamples-1
    rAbs = abs(rPos)/params.nr_samples;
    
    isl_lin = 2*sum(rAbs.^2,1)/params.nr_samples;
    isl_dB = 10*log10(isl_lin);
    
    isl_dB = permute(isl_dB, [2,3,1]);
end