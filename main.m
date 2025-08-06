close all; clearvars; clc; 
addpath("99_Utilities\utils\");

figureHandler = FigureHandler("savePath","01_Figures\");

%% Compare the mathematical equivalence between: 
% 1) Rxx(k) = sum i=1 to N-k over conj(x_i)*x_{i+k} for k =      0, ... , N-1
% 2) Rxx(k) = sum i=1 to N	 over conj(x_i)*x_{i+k} for k = -(N-1), ... , N-1

nSignal = 10;
signal = exp(1i*2*pi*randn(nSignal, 1));

[acf1, lags1] = Autocorr1(signal);
[acf2, lags2] = Autocorr2(signal);
acfIntern = xcorr(signal);
lagsIntern = lags2;

%% Plot 1
figures = {};
figures{end+1} = figure("Name","Autocorrelation1"); hold all; 
plot(lagsIntern, 20*log10(abs(acfIntern)), ...
	'-', "DisplayName","xcorr($\cdot$) function (Matlab)");

plot(lags1, 20*log10(abs(acf1)), ...
	'--', "DisplayName","$R_{xx}(k) = \sum_{n=1}^{N-k}x_n^* x_{n+k}, k\in\lbrace 0, \dots, N-1\rbrace$");

xlabel("Lags $(-)$");
ylabel("Magnitude $(dB)$");
legend("Location","northoutside", "NumColumns",2);

%% Plot 2
figures{end+1} = figure("Name","Autocorrelation2"); hold all; 
plot(lagsIntern, 20*log10(abs(acfIntern)), ...
	'-', "DisplayName","xcorr($\cdot$) function (Matlab)");

plot(lags2, 20*log10(abs(acf2)), ...
	'--', "DisplayName","$R_{xx}(k) = \sum_{n=1}^{N}x_n^* x_{n+k}, k\in\lbrace -(N-1), \dots, N-1\rbrace$");

xlabel("Lags $(-)$");
ylabel("Magnitude $(dB)$");
legend("Location","northoutside", "NumColumns",2);

%% Format Figures
for iFigure = 1:numel(figures)
	figureHandler.PreparePdf(figures{iFigure}, ...
		"width",		32,			...
		"height",		16,			...
		"font",			"Arial",	...
		"fontSize",		20,			...
		"lineWidth",	3,			...
		"markerSize",	8			...
		);
end
figureHandler.GenerateOutputs();
clear ans

%% Helper Functions
function [acf, k] = Autocorr1(signal)
	% Rxx(k) = sum i=1 to N-k over conj(x_i)*x_{i+k} for k = 0, ... , N-1
	
	N = numel(signal);
	signalPadded = paddata(signal, 2*N-1, "Side","trailing");

	acf = zeros(N, 1);
	k = 0:N-1;
	for iAcf = 1:numel(k)
		for i = 1:N-k(iAcf)
			acf(iAcf) = acf(iAcf) + conj(signalPadded(i)) * signalPadded(i+k(iAcf));
		end
	end
end

function [acf, k] = Autocorr2(signal)
	% Rxx(k) = sum i=1 to N	over conj(x_i)*x_{i+k} for k = -(N-1), ... , N-1
	
	N = numel(signal);
	signalPadded = paddata(signal, 3*N-2, "Side","both");
	% signalPadded = paddata(signalPadded, 3*N-2, "Side","trailing");

	acf = zeros(2*N-1, 1);
	k = -(N-1):(N-1);
	kBase = -k(1);
	for iAcf = 1:numel(k)
		for i = 1:N
			acf(iAcf) = acf(iAcf) + conj(signalPadded( i + kBase )) * signalPadded( i + kBase + k(iAcf));
		end
	end
end
