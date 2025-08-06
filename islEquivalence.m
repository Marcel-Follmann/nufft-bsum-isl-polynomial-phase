close all; clearvars; clc; 

addpath("99_Utilities\utils\");
figureHandler = FigureHandler("savePath","01_Figures\");

%% Numerical comparison of original ISL calculation to reformulated one

%% Parameters: 
nSamples = 20;
nParams = 3;

params = randn(nParams, 1);
paramIndexVec = (0:nParams-1).';
sampleIndexVec = 0:nSamples-1;

paramPhase = params .* (sampleIndexVec.^paramIndexVec);
signal = exp(1i*sum(paramPhase, 1));

lags = -(nSamples-1):(nSamples-1);

%% 1) ACF via xcorr
acfOriginal = xcorr(signal);

%% 2) ACF via analytical formula
% r_{xx}(k) = \sum_{n=1}^{N-k} x_n^* * x_{n+k}, k = {0, ... , N-1}

acfAnalytical = zeros(size(signal)); 
for k = 0:nSamples-1

	for n = 1:nSamples - k
		acfAnalytical(k+1) = acfAnalytical(k+1) + ...
			conj(signal(n)) * signal(n+k);
	end
end

fig = {};
fig{end+1} = figure("Name","Autocorrelation"); hold all; 

plot(lags, abs(acfOriginal), ...
	'-', "DisplayName","Original ACF");
plot(lags(nSamples:end), abs(acfAnalytical), ...
	'--', "DisplayName","Analystical ACF");

legend("Location","northeast");
