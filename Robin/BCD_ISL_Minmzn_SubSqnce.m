clc; clear; close all;
rng('shuffle');
addpath('..\common\')

%% Sequence Set Design parameters
N = 128;
L = 2;
l_optmz = 2;
Q = 3;
max_itrns = 10;
M = N/L;

lags = 1-N:N-1;

% Seed Sequence design
a_coeff_o_indx = 3;
a_coeff_i = rand(Q,L);
[xi] = create_sequence(a_coeff_i, M, L);

[ psl,isl,PSLRdB,ISLRdB ] = pslisl(xi);
InitSpecs = [ psl,isl,PSLRdB,ISLRdB ];

norm_acf_i = 10*log10(abs(xcorr(xi)).^2./N^2);
figure(1); hold all; box on;
plot(lags, norm_acf_i,'r-.');

ISL_cnvrgnce = [];
ISL_cnvrgnce(end+1) = ISLRdB;

%% Optimization loop for chosen a_coeff [i.e. first, second or third a_coeff]
up_lim = 100;
lwr_lim = -100;
a_points = 1000;
a_var = linspace(lwr_lim, up_lim, a_points);

a_coeff = a_coeff_i;
for itr = 1: max_itrns
for q=1:Q
a_coeff_o_indx = q;
ISL_a_cntnr = zeros(2,length(a_var));

for i = 1:length(a_var)
    a_coeff(a_coeff_o_indx, l_optmz) = a_var(i);
    [xv] = create_sequence(a_coeff, M, L);

    [ psl,isl,PSLRdB,ISLRdB ] = pslisl(xv);

    ISL_a_cntnr(:,i) = [isl; a_var(i)];
end
[~,indx] = min(ISL_a_cntnr(1,:));

a_coeff(a_coeff_o_indx, l_optmz) = ISL_a_cntnr(2, indx);
end

%% Capture the ISL results for current iteration
[x_int] = create_sequence(a_coeff, M, L);
[ psl,isl,PSLRdB,ISLRdB ] = pslisl(x_int);
ISL_cnvrgnce(end+1) = ISLRdB;
end


a_coeff_fnl = a_coeff;
xf = x_int;
[ psl,isl,PSLRdB,ISLRdB ] = pslisl(xf);
Final_Specs = [ psl,isl,PSLRdB,ISLRdB ];
norm_acf_f = 10*log10(abs(xcorr(xf)).^2./N^2);
figure(1);
plot(lags, norm_acf_f,'b-');


figure(2);
plot(ISL_cnvrgnce);
xlabel('iterations')
ylabel('obj')

%% Support functions
function [x] = create_sequence(a_coeff, M, L)
theta_li = [];
m = (1:M).';
for l = 1:L
    theta_li(:,l) = a_coeff(1,l)+ a_coeff(2,l).*m + a_coeff(3,l).*m.^2;
end
theta_x = theta_li(:);
x = exp(1j*theta_x);
end