close all; clearvars; clc; 

N = 20;
n = 0:N-1;
A = pi;
da = pi/N;
os = 100;
a = (0:da/os:A-da/os).';

x = exp(1i*(a.*n.^2));

r_xx = zeros(N, 2*N-1);
for k = 1:numel(a)
	r_xx(k,:) = xcorr(x(k,:));
end

isl = sum(abs(r_xx).^2,2);



plot(a/pi,20*log10(isl));
xticks((0:da:A-da)/pi);
