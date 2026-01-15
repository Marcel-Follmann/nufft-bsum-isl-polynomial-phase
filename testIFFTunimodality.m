close all; clearvars; clc; 

addpath("01_Figures\");
addpath("99_Utilities\utils\");
figureHandler = FigureHandler("savePath","01_Figures\");

%% Simulate a frequency spectrum -> convert it to time -> evaluate if there is a unimodality guarantee between to adjacent samples

t = linspace(0,0.5);
F = [1, 1/3];
x = cos(2*pi*F(1)*t) + 0.5*cos(2*pi*F(2)*t);

plot(t,x)