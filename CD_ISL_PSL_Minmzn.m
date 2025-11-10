close all; clearvars; clc; 
rng('shuffle');

addpath("99_Utilities\utils\");

figureHandler = FigureHandler("savePath","01_Figures\");

% rng(1);

%% Sequence Set Design parameters
nSequence = 20;
paramOrder = 4;
nIterations = 5;

lags = 1-nSequence:nSequence-1;

% flag for chosing objective of interest (f_ooi) - option1: 'PSL'/ option2: 'ISL'
objectiveFlag = 'ISL';   % 'PSL' / 'ISL'

objectiveConvegence = [];

% Seed Sequence design
paramInitial = [-0.920274947429725;1.179610379138849;1.508569156483106;1.062646144475084]; %rand(paramOrder,1);
[sequenceInitial] = CreateSequence(paramInitial, nSequence, paramOrder);

[ psl,isl,pslr_dB,islr_dB ] = pslisl(sequenceInitial);
InitSpecs = [ psl,isl,pslr_dB,islr_dB ];
switch objectiveFlag
    case 'PSL'
        objectiveConvegence(1) = pslr_dB;
    case 'ISL'
        objectiveConvegence(1) = islr_dB;
end



%% Optimization loop for chosen a_coeff [i.e. first, second or third a_coeff]
paramLimitUpper = 1000;
paramLimitLower = -1000;
nParamVec = 1000;
paramVec = linspace(paramLimitLower, paramLimitUpper, nParamVec);

paramOptimized = paramInitial;

% iterate over multiple iterations to optimize the objective
for iIteration = 1: nIterations
	% iterate over all the phase coefficients
	for iOrder=1:paramOrder
		objectiveParamArray = zeros(2,length(paramVec));
		% iterate over the entire range of phase coefficient {lwr_lim, up_lim}
		for iParam = 1:length(paramVec)
    		paramOptimized(iOrder) = paramVec(iParam);
    		[sequenceOptimized] = CreateSequence(paramOptimized, nSequence, paramOrder);
		
    		[ psl,isl,~, ~ ] = pslisl(sequenceOptimized);
		
    		switch objectiveFlag
        		case 'PSL'
            		objectiveParamArray(:,iParam) = [psl; paramVec(iParam)];
        		case 'ISL'
            		objectiveParamArray(:,iParam) = [isl; paramVec(iParam)];
    		end
		end
		[~,indx] = min(objectiveParamArray(1,:));
		
		paramOptimized(iOrder) = objectiveParamArray(2, indx);
	end

%% Capture the ISL results for current iteration
[sequenceOptimized] = CreateSequence(paramOptimized, nSequence, paramOrder);
[ ~,~,pslr_dB,islr_dB ] = pslisl(sequenceOptimized);
	switch objectiveFlag
    	case 'PSL'
        	objectiveConvegence(iIteration + 1) = pslr_dB;
    	case 'ISL'
        	objectiveConvegence(iIteration + 1) = islr_dB;
	end
end


%% Plot
fig = {};
[ psl,isl,pslr_dB,islr_dB ] = pslisl(sequenceOptimized);
Final_Specs = [ psl,isl,pslr_dB,islr_dB ];

acfOptimized = 10*log10(abs(xcorr(sequenceOptimized)).^2./nSequence^2);
acfInitial = 10*log10(abs(xcorr(sequenceInitial)).^2./nSequence^2);

fig{end+1} = figure("Name", "Autocorrelation Order " + num2str(paramOrder)); hold all; 
plot(lags, acfOptimized,'b-', "DisplayName","Optimized Sequence");
plot(lags, acfInitial,'r-.', "DisplayName","Initial Sequence");
xlabel('lag $(-)$')
ylabel('auto-correlation $(-)$')
legend("Location","best");


fig{end+1} = figure("Name","Objective Order" + num2str(paramOrder));
plot(0:nIterations,objectiveConvegence);
xlabel('iterations')
ylabel('obj')

phaseInitial = unwrap(atan2(imag(sequenceInitial),real(sequenceInitial)));
phaseOptimized = unwrap(atan2(imag(sequenceOptimized),real(sequenceOptimized)));

fig{end+1} = figure("Name", "Phases Order" + num2str(paramOrder)); hold all; 
indx = (1:nSequence)';
plot(indx, phaseInitial, '-.ro', "DisplayName","Initial Phase");
plot(indx, phaseOptimized, '-.bo', "DisplayName","Optimized Phase");
xlabel('Code length');
ylabel('Phase');
legend("Location","best");


for iFig = 1:numel(fig)
	figureHandler.PreparePdf(fig{iFig}, ...
		"width",		32,			...
		"height",		16,			...
		"font",			"Arial",	...
		"fontSize",		20,			...
		"lineWidth",	3,			...
		"markerSize",	8			...
		);

end
figureHandler.GenerateOutputs();

disp(paramOptimized)

%% Support functions
function [sequence] = CreateSequence(params, nSequence, paramOrder)
    iSequenceVec = (1:nSequence).';
    % phase = params(1) + params(2).*iSequenceVec + params(3).*iSequenceVec.^2;
    phase = 0;
    for iOrder=1:paramOrder
        phase = phase + params(iOrder).*iSequenceVec.^(iOrder-1);
	end

    sequence = exp(1j*phase(:));
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
isl         = 2*sum(r(1:N-1).^2);

PSLRdB      = 10*log10(psl^2/N^2);
ISLRdB      = 10*log10(isl/N^2);
end

