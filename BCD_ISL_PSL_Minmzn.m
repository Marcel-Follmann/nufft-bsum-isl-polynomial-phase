clc; clear; close all;
rng('shuffle');
% rng(1);
% addpath('..\common\')
%% Sequence Set Design parameters
N = 128;
Q = 3;
max_itrns = 5;

lags = 1-N:N-1;

% flag for chosing objective of interest (f_ooi) - option1: 'PSL'/ option2: 'ISL'
f_ooi = 'ISL';   % 'PSL' / 'ISL'

Obj_cnvrgnce = [];

% Seed Sequence design
a_coeff_o_indx = 3;
a_coeff_i = rand(Q,1);
[xi] = create_sequence(a_coeff_i, N, Q);

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

% iterate over multiple iterations to optimize the objective
for itr = 1: max_itrns
% iterate over all the phase coefficients
for q=1:Q
a_coeff_o_indx = q;
Obj_a_cntnr = zeros(2,length(a_var));
% iterate over the entire range of phase coefficient {lwr_lim, up_lim}
for i = 1:length(a_var)
    a_coeff(a_coeff_o_indx) = a_var(i);
    [xv] = create_sequence(a_coeff, N, Q);

    [ psl,isl,PSLRdB,ISLRdB ] = pslisl(xv);

    switch f_ooi
        case 'PSL'
            Obj_a_cntnr(:,i) = [psl; a_var(i)];
        case 'ISL'
            Obj_a_cntnr(:,i) = [isl; a_var(i)];
    end
end
[~,indx] = min(Obj_a_cntnr(1,:));

a_coeff(a_coeff_o_indx) = Obj_a_cntnr(2, indx);
end

%% Capture the ISL results for current iteration
[x_int] = create_sequence(a_coeff, N, Q);
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
xlabel('lag')
ylabel('auto-correlation')
legend('init seq','opt seq')


figure(2);
plot(Obj_cnvrgnce);
xlabel('iterations')
ylabel('obj')

arg_x = unwrap(atan2(imag(xi),real(xi)));
arg_xf = unwrap(atan2(imag(xf),real(xf)));
figure(24); box on; hold all; % grid on;
indx = (1:N)';
plot(indx, arg_x, '-.ro', 'linewidth', 1.5);
plot(indx, arg_xf, '-.bo', 'linewidth', 1.5);
xlabel('Code length');ylabel('Phase');
hold off;

disp(a_coeff_fnl)
%% Support functions
function [x] = create_sequence(a_coeff, M, Q)
    m = (1:M).';
    theta_l = a_coeff(1) + a_coeff(2).*m + a_coeff(3).*m.^2;
    theta_li = 0;
    for q=1:Q
        theta_li = theta_li + a_coeff(q).*m.^(q-1);
    end
    theta_x = theta_li(:);

    x = exp(1j*theta_x);
end



function [ psl,isl,PSLRdB,ISLRdB ] = pslisl( x , nF , pF )
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Copyright (c) 2020, <Mohammad Alaee>, website: https://radarmimo.com/
% All rights reserved.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 3. All advertising materials mentioning features or use of this software
%    must display the following acknowledgement:
%    This product includes software developed by the <organization>.
% 4. Neither the name of the <organization> nor the
%    names of its contributors may be used to endorse or promote products
%    derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY <COPYRIGHT HOLDER> ''AS IS'' AND ANY
% EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Website: https://radarmimo.com/
% Email: info@radarmimo.com, mohammad.alaee@uni.lu, mmalaee@yahoo.com
% Code written by : Mohammad Alaee
% Last update : 2016
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This function calculates PSL and ISL for a radar transmit sequence 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Input: 
%           x       : sequence
%           nF      : normalization flag, default n = 0 
%           pF      : pF = 1 -> calculating "Periodic" PSL and ISL,
%                     otherwise, calculating "Aperiodic" PSL and ISL,
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Output: 
%           psl     : calculated PSL value
%           isl     : calculated ISL value
%           PSLRdB  : calculated PSLR (dB) value
%           ISLRdB  : calculated ISLR (dB) value
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% The input parameters can be chosed as follows:
if nargin == 1
    nF      = 0;
    pF      = 0; 
elseif nargin <= 2
   pF       = 0; 
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x           = x(:);
N           = length(x);   
if pF == 0
    r       = abs(xcorr(x));
else
    r       = abs(pcorr(x));
end
if nF == 1
    r       = r / max(r);
end
psl         = max(r(1:N-1));
isl         = sum(r(1:N-1).^2);

PSLRdB      = 10*log10(psl^2/N^2);
ISLRdB      = 10*log10(isl/N^2);
end

