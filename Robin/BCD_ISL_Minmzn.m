clc; clear; close all;
% addpath('C:\Program Files\Mosek\10.2\toolbox\r2017aom')
%% Sequence Set Design parameters
N = 128;
L = 1;
M = N/L;
rng(1);
addpath('..\common\')
epsln = 1e-3;  % epsln = 1e-3;
scale_fctr = 1;
Q = 3;
max_itrns = 30;
theta_li = [];
m = (1:M).';
a_coeff = zeros(Q,L);
for l = 1:L
    a_coeff(:,l) = scale_fctr*rand(Q,1);
    theta_li(:,l) = a_coeff(1,l) + a_coeff(2,l).*m + a_coeff(3,l).*m.^2;
end
theta_x = theta_li(:);
x = exp(1j*theta_x);

[ psl,isl,PSLRdB,ISLRdB ] = pslisl(x)

norm_acf = 10*log10(abs(xcorr(x)).^2./N^2);

figure;
plot(norm_acf);


up_lim = 10;
lwr_lim = -10;
a_points = 1000;
a_var = linspace(lwr_lim, up_lim, a_points);

l = 1;
for i = 1:length(a_var)
    a_coeff(1,l) = a_var(i);
    theta_li(:,l) = a_coeff(1,l) + a_coeff(2,l).*m + a_coeff(3,l).*m.^2;
    theta_x = theta_li(:);
end