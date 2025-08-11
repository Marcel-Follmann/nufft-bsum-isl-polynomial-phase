close all; clearvars; clc;
load("g.mat")

x = ifftshift(g);
X = ifft(x);
X = fftshift(X);
plot(X)