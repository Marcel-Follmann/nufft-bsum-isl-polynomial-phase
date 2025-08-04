close all; clearvars; clc; 
addpath("99_Utilities\utils\");

rng default

a = ones(3,1);

N = 10;
n = 0:N-1;

a(3) = pi/N+pi;

x = exp(1i*(a(1)+a(2)*n + a(3)*n.^2));
r_xx = xcorr(x);
isl = sum(abs(r_xx).^2);
