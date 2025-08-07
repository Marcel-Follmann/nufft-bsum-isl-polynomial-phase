close all; clearvars; clc; 

addpath("99_Utilities\utils\");
figureHandler = FigureHandler("savePath","01_Figures\");

%% Numerical comparison of original ISL calculation to reformulated one

%% Parameters: 
nSamples = 5;
nParams = 3;

d = zeros(nSamples-1, nSamples - 1, nSamples);
for k = 0:nSamples-1
	for n = 1:nSamples - k
		for m = 1:nSamples - k
			d(n,m,k+1) = eqn.deltaD(n,m,k,5);
		end
	end
	% surf(d(:,:,k+1));
	% zlim(200*[-1,1])
	% pause(.5);
end

% Flatten your 3D array and remove the zeros
D = d(:);
D = D(D~=0);

% Count occurrences of each distinct d-value
[vals, ~, idx] = unique(D);
counts = accumarray(idx,1);

% Now plot with plot()
figure;
plot(vals, counts, '-o','LineWidth',1.5);
xlabel('d value');
ylabel('Number of occurrences');
grid on;



