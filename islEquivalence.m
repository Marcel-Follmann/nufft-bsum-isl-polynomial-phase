close all; clearvars; clc; 

addpath("99_Utilities\utils\");
figureHandler = FigureHandler("savePath","01_Figures\");

%% Numerical comparison of original ISL calculation to reformulated one

%% Parameters: 
numSamples = 20;
numParams = 3;
p = numParams;

params = randn(numParams, 1);
paramIndexVec = (0:numParams-1).';
sampleIndexVec = 0:numSamples-1;

phaseParam = params .* (sampleIndexVec.^paramIndexVec);
signal = exp(1i*sum(phaseParam, 1));

lags = -(numSamples-1):(numSamples-1);

%% 1) ISL via xcorr
acfOriginal = xcorr(signal);% / numSamples;
islOriginal = Isl(acfOriginal);

%% 2) ISL via analytical double-sum
islAnalytic = 0;
for k = 1:(numSamples-1)
	
	% Constant term
	termConstant = 0;
	for n = 1:numSamples
		termConstant = termConstant + abs(ComputeResidualComplexWeight(n,k,params,p)).^2;
	end

	% Cosine term
	termCosine = 0;
	for n=1:numSamples
		for m=n+1:numSamples
			termCosine = termCosine + ...
				ComputeAmplitude(n,m,k,p,params) * ...
				cos( ...
					ComputeFrequencyIndex(n,m,k,p) * params(p) + ...
					ComputePhaseOffset(n,m,k,params,p) ...
				);
		end
	end

	islAnalytic = islAnalytic + 2*(termConstant + termCosine);
end

fig = {};
fig{end+1} = figure("Name","Autocorrelation"); hold all; 

plot(lags, 20*log10(abs(acfOriginal)));


% Original ISL calculation


function isl_linear = Isl(autocorrelation)
	numSamples = (numel(autocorrelation)+1) / 2;

	isl_linear = sum(abs(autocorrelation).^2) - abs(autocorrelation(numSamples));
	isl_dB = 10*log10(isl_linear);
end

