close all; clearvars; clc; 

addpath("..\99_Utilities\utils\");
figureHandler = FigureHandler("savePath","..\01_Figures\");
fig = {};

%% Equal Spaced Frequency -> Non-equal Spaced Time
N = 50;
K = 100;
J = 20;
M = 1000;

s = ones(N,1);

nufft = NuFft("J",J, "s", s);

omega_m = 2*pi*sort(rand(1,M)-0.5);

% Frequency Domain
A = 1;
n(:,1) = -N/2:N/2-1; 
xn(:,1) = A/4*(n == 10) + A/4*(n == -10);

% Transform
Xn(:,1) = fft(fftshift(xn));
fn(:,1) = (-1/2:1/N:1/2-1/N);

[Y, f] = nufft.kPointFft(xn, K);
f = f;
% Interpolation
nufft.prepareInterpolator();
for m = 1:numel(omega_m)
    u = nufft.interpolators(omega_m(m));
    j = 1:J;
    idx = mod( nufft.k0( omega_m(m) ) + j , K ) + 1;
    Xm(m) = sum( Y( idx ) .* conj(u) );
end



%% Visualization
fig{end+1} = figure("Name","Frequency-Domain"); hold on; 
stem(n, xn);
box on; 

fig{end+1} = figure("Name","Time-Domain"); hold on; 
plot(fn, real(Xn), 'o', 'DisplayName','N-Point FFT');
plot(f, real(Y), '.', 'DisplayName', 'K-Point FFT');
plot(omega_m/(2*pi), real(Xm), 'o', 'DisplayName', 'NUFFT');

box on; 
legend("Location","northeast")


for i = 1:numel(fig)
    figureHandler.PrepareEmf(fig{i});
end
