close all; clearvars; clc; 

n = 100;
f0 = 10;
fOs = 2;
nOs = fOs*n;

X = zeros(n,1);
f = -n/2 : n/2 - 1;
idx0 = n/2+1;
X(idx0 + f0*[1,-1]) = n/2;

fMod = f; 
XMod = X;
fMod(XMod == 0) = [];
XMod(XMod == 0) = [];

t = 0:1/(2*n):1-1/(2*n);
tOs = 0:1/nOs:1-1/nOs;
x = ifft(ifftshift(X), nOs);
xMod = nufft(XMod, fMod, t)/nOs; 

nexttile; hold on; 
plot(f, X);
plot(fMod, XMod, 'o');
xlabel("Frequency (Hz)");

nexttile; hold on; 
plot(tOs, x)
plot(t, xMod, 'o');
xlabel("Time (s)");