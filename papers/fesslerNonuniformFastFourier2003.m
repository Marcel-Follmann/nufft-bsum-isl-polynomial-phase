close all; clearvars; clc; 

addpath("..\99_Utilities\utils\");
figureHandler = FigureHandler("savePath","..\01_Figures\");
fig = {};

%% Equal Spaced Time -> Nonequal Spaced Frequency
N = 100;
K = 1000;
A = 1;
n(:,1) = -N/2:N/2-1; 
xn(:,1) = A/4*(n == 10) + A/4*(n == -10);

Xn(:,1) = fft(fftshift(xn));
for k = 0:K/N-1
    Xnk(:,k+1) = fft(fftshift(xn.*exp(-1i*2*pi/K*n*k)));
end
Xnk = reshape(Xnk.', [], 1);
f(:,1) = -1/2:1/N:1/2-1/N;
fk(:,1) = -1/2:1/K:1/2-1/K; 

fig{end+1} = figure("Name","Frequency-Domain"); hold on; 
stem(n, xn);
box on; 

fig{end+1} = figure("Name","Time-Domain"); hold on; 
% plot(f, Xn, '.');
plot(f, real(Xn), 'o', 'DisplayName','N-Point FFT');
plot(fk, real(Xnk), '.', 'DisplayName', 'K-Point FFT');

box on; 
legend("Location","northeast")


for i = 1:numel(fig)
    figureHandler.PrepareEmf(fig{i});
end
