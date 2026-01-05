close all; clearvars; clc; 

addpath("..\99_Utilities\utils\");
figureHandler = FigureHandler("savePath","..\01_Figures\");
fig = {};

%% Parameters

k(:,1) = 10*[-1,1];
f(1:2, 1) = .5; 

X = 40;

dx = .125; 
ds = 1; 
Mtau = 2*pi/(dx*ds);

%% Simulation
g = @(x, tau) exp(-x.^2/(4*tau));

ndx = -X/2:dx:X/2-dx;
fTau = sum(f.*g(ndx - k, dx*10), 1);

%% Visualization
fig{end+1} = figure("Name", "Frequency Domain"); hold on; 
stem(k, f);
plot(ndx, fTau, '.');
