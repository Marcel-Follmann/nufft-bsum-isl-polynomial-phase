close all; clearvars; clc; 


f0 = 20;
f1 = 40;
fs = 200;
fOs = 200;

n = fs; 

nOs = fOs;

X = zeros(n,1);
f = -n/2 : n/2 - 1;
idx0 = n/2+1;
X(idx0 + f0*[1,-1]) = exp(1i*randn(1)*[-1,1]) * n/2;
X(idx0 + f1*[1,-1]) = exp(1i*randn(1)*[-1,1]) * n/2;

fMod = f; 
XMod = X;
fMod(XMod == 0) = [];
XMod(XMod == 0) = [];

t = 0:1/(n):1-1/(n);
tMod = 0:1/nOs:1-1/nOs;
x = ifft(ifftshift(X));
xMod = nufft(flip(XMod), fMod, tMod)/n; 

nexttile; hold on; 
plot(f, X);
plot(fMod, XMod, 'o');
xlabel("Frequency (Hz)");

nexttile; hold on; 
plot(t, x)
plot(tMod, xMod, '-o');
xlabel("Time (s)");