classdef TestMathematicalDerivations < matlab.unittest.TestCase
	properties
		nSamples
		nParams
		params
		signal
		lags
		acfOriginal
	end

	methods (TestClassSetup)
		function SetupClass(tc)
			close all; clearvars; clc; 

			tc.nSamples = 20;
			tc.nParams = 3;
			
			tc.params = randn(tc.nParams, 1);
			paramIndexVec = (0:tc.nParams-1).';
			sampleIndexVec = 0:tc.nSamples-1;
			
			paramPhase = tc.params .* (sampleIndexVec.^paramIndexVec);
			tc.signal = exp(1i*sum(paramPhase, 1));
			
			tc.lags = 0:(tc.nSamples-1);
			
			%% 1) ACF via xcorr
			acfOriginal_full = xcorr(tc.signal);
			tc.acfOriginal = acfOriginal_full(tc.nSamples:end);
		end
	end

	methods (TestMethodSetup)
		% Setup for each test
	end
		
	methods (Test)
		% Test methods

		function Test_Autocorrelation(tc) 
			acfAnalytical = zeros(size(tc.acfOriginal)); 
			for k = 0:tc.nSamples-1
			
				for n = 1:tc.nSamples - k
					acfAnalytical(k+1) = acfAnalytical(k+1) + ...
						conj(tc.signal(n)) * tc.signal(n+k);
				end
			end
			plot(tc.lags, abs(tc.acfOriginal)); hold on;
			% plot(tc.lags, abs(acfAnalytical));
			pause(10);
		end
	end

end