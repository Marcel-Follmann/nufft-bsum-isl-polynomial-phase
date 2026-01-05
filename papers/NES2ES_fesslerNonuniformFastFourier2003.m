close all; clearvars; clc; 

addpath("..\99_Utilities\utils\");
figureHandler = FigureHandler("savePath","..\01_Figures\");
fig = {};

%% Equal Spaced Frequency -> Non-equal Spaced Time
sigma =2;
N = 1000; 
w = 50;
alpha = pi*(sigma-1) / w;

n(1,:) = (-sigma*N/2:sigma*N/2 - 1);
M = numel(n);

% Frequency Domain
A = 1;
x(:,1) = 10*[-1,1];
f(:,1) = A/4 * ones(2,1);

phi = @(x, alpha) exp(-alpha*x.^2).* (abs(x) <= w);        % kernel space-domain
Phi = @(x, alpha) sqrt(pi/alpha) * exp(-(2*pi*x).^2/(4*alpha));    % kernel frequency-domain

g = sum(f.*phi(x - n, alpha), 1); % Convolution with kernel
G = fftshift(fft(ifftshift(g), M)); % Transform

k = -M/2:(M/2-1);
omega = k/M;


% keep central N bins
iKeep = floor((M-N)/2);
G = G(iKeep+1 : iKeep+N);
omega = omega(iKeep+1 : iKeep+N);

F = G./Phi(omega, alpha);



%% Visualization
fig{end+1} = figure("Name","Frequency-Domain"); hold on; 
plot((x* [1, 1]).', (f*[0, 1]).') 
plot(n, g, '-o');
xlim(2*x)
box on; 
xlabel("Frequency (Hz)");
ylabel("Magnitude (-)");

fig{end+1} = figure("Name","Time-Domain"); hold on; 
plot(omega, real(F), '-o');
xlabel("Time (s)");
ylabel("Amplitude (-)");
% xlim([-1, 1])

for i = 1:numel(fig)
    figureHandler.PrepareEmf(fig{i});
end
