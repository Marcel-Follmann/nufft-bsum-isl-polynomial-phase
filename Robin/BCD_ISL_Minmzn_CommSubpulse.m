clc; clear; close all;
rng('shuffle');
% rng(1);
addpath('..\common\')

%% Sequence Set Design parameters
N = 500;
L = 2;
l_optmz = 2;    % Sensing subpulse index for optimization
l_comm = 1;     % Communication subpulse index
Q = 3;
max_itrns = 10;
M = N/L;
% flag for chosing objective of interest (f_ooi) - option1: 'PSL'/ option2: 'ISL'
f_ooi = 'ISL';  
lags = 1-N:N-1;
symbl = 3;      % Symbol input range: 0 to 7

% Seed Sequence design
a_coeff_o_indx = 3;
a_coeff_i = rand(Q,L);

[a_coeff_i] = comm_subpulse_coeff(symbl, a_coeff_i, l_comm);

[xi] = create_sequence(a_coeff_i, M, L);

[ psl,isl,PSLRdB,ISLRdB ] = pslisl(xi);
InitSpecs = [ psl,isl,PSLRdB,ISLRdB ];

norm_acf_i = 10*log10(abs(xcorr(xi)).^2./N^2);
figure(1); hold all; box on;
plot(lags, norm_acf_i,'r-.');

Obj_cnvrgnce = [];
switch f_ooi
    case 'PSL'
        Obj_cnvrgnce(end+1) = PSLRdB;
    case 'ISL'
        Obj_cnvrgnce(end+1) = ISLRdB;
end


%% Optimization loop for chosen a_coeff [i.e. first, second or third a_coeff]
up_lim = 100;
lwr_lim = -100;
a_points = 1000;
a_var = linspace(lwr_lim, up_lim, a_points);

a_coeff = a_coeff_i;
for itr = 1: max_itrns
for q=1:Q
a_coeff_o_indx = q;
Obj_a_cntnr = zeros(2,length(a_var));

for i = 1:length(a_var)
    a_coeff(a_coeff_o_indx, l_optmz) = a_var(i);
    [xv] = create_sequence(a_coeff, M, L);

    [ psl,isl,~,~ ] = pslisl(xv);
    
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

a_coeff = a_coeff_fnl;

[ISL_comm] = compute_ISL_Comm(a_coeff, Q, symbl, l_comm, M, L);


%% Support functions

% Generate a sequence from phase coefficients
function [x] = create_sequence(a_coeff, M, L)
theta_li = [];
m = (1:M).';
for l = 1:L
    theta_li(:,l) = a_coeff(1,l)+ a_coeff(2,l).*m + a_coeff(3,l).*m.^2;
end
theta_x = theta_li(:);
x = exp(1j*theta_x);
end

% Assign phase coefficients depending upon the chosen symbol for communication 
function [a] = comm_subpulse_coeff(symbl, a, lc)
switch symbl
    case 0
        a(:,lc) = [-1, -1, -1];
    case 1
        a(:,lc) = [-1,  1, -1];
    case 2
        a(:,lc) = [ 1, -1, -1];
    case 3
        a(:,lc) = [ 1,  1, -1];
    case 4
        a(:,lc) = [-1, -1,  1];
    case 5
        a(:,lc) = [-1,  1,  1];
    case 6
        a(:,lc) = [ 1, -1,  1];
    case 7
        a(:,lc) = [ 1,  1,  1];
end
end

% sc -> chosen symbol
function [ISL_comm] = compute_ISL_Comm(a, Q, sc, l_comm, M, L)
num_symb = 2^Q;
ISL_comm = zeros(num_symb,1);
for j = 1:num_symb
    symb = j-1;
%     if j-1 == sc
%         continue;
%     end
    [a] = comm_subpulse_coeff(symb, a, l_comm);
    [x] = create_sequence(a, M, L);
    [ ~,~,~,ISL_comm(j) ] = pslisl(x);
end
end
