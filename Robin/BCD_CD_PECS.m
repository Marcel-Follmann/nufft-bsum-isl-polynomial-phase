clc; clear; close all;
rng('shuffle');
% addpath('..\common\')

%% Sequence Set Design parameters
N = 64;
L = 1;

Q = 3;
max_itrns = 5;
M = N/L;

lags = 1-N:N-1;
% flag for chosing objective of interest (f_ooi) - option1: 'PSL'/ option2: 'ISL'
f_ooi = 'ISL';  

% Seed Sequence design
a_coeff_o_indx = 3;
a_coeff_i = rand(Q,L);
[xi] = create_sequence(a_coeff_i, M, L);

Obj_cnvrgnce = [];

[ psl,isl,PSLRdB,ISLRdB ] = pslisl(xi);
InitSpecs = [ psl,isl,PSLRdB,ISLRdB ];
switch f_ooi
    case 'PSL'
        Obj_cnvrgnce(1) = PSLRdB;
    case 'ISL'
        Obj_cnvrgnce(1) = ISLRdB;
end

norm_acf_i = 10*log10(abs(xcorr(xi)).^2./N^2);
figure(1); hold all; box on;
plot(lags, norm_acf_i,'r-.');




%% Optimization loop for chosen a_coeff [i.e. first, second or third a_coeff]
up_lim = 1000;
lwr_lim = -1000;
a_points = 10000;
a_var = linspace(lwr_lim, up_lim, a_points);

a_coeff = a_coeff_i;
for itr = 1: max_itrns
    for l = 1:L
        l_optmz = l;
        for q=1:Q
            a_coeff_o_indx = q;
            Obj_a_cntnr = zeros(2,length(a_var));

            for i = 1:length(a_var)
                a_coeff(a_coeff_o_indx, l_optmz) = a_var(i);
                [xv] = create_sequence(a_coeff, M, L);

                [ psl,isl,PSLRdB,ISLRdB ] = pslisl(xv);

                switch f_ooi
                    case 'PSL'
                        Obj_a_cntnr(:,i) = [psl; a_var(i)];
                    case 'ISL'
                        Obj_a_cntnr(:,i) = [isl; a_var(i)];
                end
            end
            [~,indx] = min(Obj_a_cntnr(1,:));

            a_coeff(a_coeff_o_indx, l_optmz) = Obj_a_cntnr(2, indx);
        end

    end
%% Capture the ISL results for current iteration
[x_int] = create_sequence(a_coeff, M, L);
[ psl,isl,PSLRdB,ISLRdB ] = pslisl(x_int);
switch f_ooi
    case 'PSL'
        Obj_cnvrgnce(end+1) = PSLRdB;
    case 'ISL'
        Obj_cnvrgnce(end+1) = ISLRdB;
end
end


a_coeff_fnl = a_coeff;
xf = x_int;
[ psl,isl,PSLRdB,ISLRdB ] = pslisl(xf);
Final_Specs = [ psl,isl,PSLRdB,ISLRdB ];
norm_acf_f = 10*log10(abs(xcorr(xf)).^2./N^2);
figure(1);
plot(lags, norm_acf_f,'b-');


figure(2);
plot(Obj_cnvrgnce);
xlabel('iterations')
ylabel('obj')



figure; hold all; box on;


indx_e = 0;
% iterate over subsequences
for l = 1:L
    % -- indexing --
    indx_s = indx_e+1;
    indx_e = indx_s + M-1;
    cp = xf(indx_s:indx_e);
    psi_cp = unwrap(atan2(imag(cp),real(cp)));
    plot(indx_s:indx_e, psi_cp, '--o')
    xlabel('code index')
    ylabel('phase [rad]')
end
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