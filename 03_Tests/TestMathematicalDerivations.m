classdef TestMathematicalDerivations < matlab.unittest.TestCase
	properties
		nSamples
		nParams
		params
		signal
		lags
		acfOriginal
		islOriginal
		relTol = 1e-3
	end

	methods (TestClassSetup)
		function SetupClass(tc)
			close all; 
			
			tc.nSamples = 20;
			tc.nParams = 4;
			
			tc.params = randn(tc.nParams, 1);
			paramIndexVec = (0:tc.nParams-1).';
			sampleIndexVec = 0:tc.nSamples-1;
			
			paramPhase = tc.params .* (sampleIndexVec.^paramIndexVec);
			tc.signal = exp(1i*sum(paramPhase, 1));
			
			tc.lags = 0:(tc.nSamples-1);
			
			acfOriginal_full = xcorr(tc.signal);
			tc.acfOriginal = acfOriginal_full(tc.nSamples:end);
			tc.islOriginal = 2*sum( abs(tc.acfOriginal(2:end)).^2 );
		end
	end

	methods (TestMethodSetup)
		
	end
		
	methods (Test)
		% Test methods

		function Test_AutocorrelationDefinition(tc) 
			% r_xx(k) = \sum_{n=1}^{N-k} x_n^* * x_{n+k}, k = {0, 1, ..., N-1}
			
			acfAnalytical = zeros(size(tc.acfOriginal)); 
			for k = 0:tc.nSamples-1
			
				for n = 1:tc.nSamples-k
					acfAnalytical(k+1) = acfAnalytical(k+1) + ...
						conj(tc.signal(n)) * tc.signal(n+k);
				end
			end
			
			% tc.PlotAcf(acfAnalytical);

			tc.verifyEqual(abs(acfAnalytical), abs(tc.acfOriginal), "RelTol", tc.relTol);
		end

		function Test_AutocorrelationSimplifiedComplete(tc) 
			% r_xx(k; a) = \sum_{n=1}^{N-k} exp{ j [ \sum_{q=1}^{Q} a_q * d_{n,k}^{(q)} ] }

			acfAnalytical = zeros(size(tc.acfOriginal));
			for k = 0:tc.nSamples-1
				
				for n = 1:tc.nSamples-k
					% \sum_{q=1}^{Q} a_q * d_{n,k}^{(q)}
					phaseSum = 0;
					for q=1:tc.nParams-1 % a_0, a_1, ..., a_{Q-1}
						phaseSum = phaseSum + ...
							tc.params(q+1) * eqn.d(n,k,q);
					end

					acfAnalytical(k+1) = acfAnalytical(k+1) + ...
						exp(1j*phaseSum);
				end
			end

			tc.PlotAcf(acfAnalytical);

			tc.verifyEqual(abs(acfAnalytical), abs(tc.acfOriginal), "RelTol", tc.relTol);
		end

		function Test_AutocorrelationSimplifiedProblemSplit(tc) 
			% r_xx(k; a_p, a) = \sum_{n=1}^{N-k} exp{ j a_p * d_{n,k}^{(p)} } * c_{n,k}^{(p)}
			
			p = 2; %randi(tc.nParams-1); % select a_1 to a_{Q-1}

			acfAnalytical = zeros(size(tc.acfOriginal));
			for k = 0:tc.nSamples-1
				
				for n = 1:tc.nSamples-k
					acfAnalytical(k+1) = acfAnalytical(k+1) + ...
						exp( 1j* tc.params(p) * eqn.d(n,k,p) ) * eqn.c(n,k,tc.params, p);
				end
			end

			% tc.PlotAcf(acfAnalytical);

			tc.verifyEqual(abs(acfAnalytical), abs(tc.acfOriginal), ...
				sprintf("Failed for p = %d", p), ...
				"RelTol", tc.relTol);
		end

		function Test_IslProblemSplit(tc)
			% \Phi_{ISL}(a_p;a) = 2\sum_{k=1}^{N-1} \sum_{n=1}^{N-k}
			%		\sum_{m=1}^{N-k} c_{n,k}^{(p)}(c_{m,k}^{(p)})^* 
			%		* exp(j a_p \delta d_{n,m,k}^{(p)})

			p = randi(tc.nParams-1); % select a_1 to a_{Q-1}

			islAnalytical = 0;
			for k = 1:tc.nSamples-1
				
				for n=1:tc.nSamples-k
					for m=1:tc.nSamples-k
						islAnalytical = islAnalytical + ...
							2*eqn.c(n,k,tc.params,p) * ...
							conj(eqn.c(m,k,tc.params,p)) * ...
							exp(1j*tc.params(p) * eqn.deltaD(n,m,k,p));
					end
				end
				
			end

			tc.verifyEqual(islAnalytical, tc.islOriginal);

		end

	end

	methods
		function PlotAcf(tc, acf)
			figure("Name","Autocorrelation"); hold all; 
			plot(tc.lags, abs(tc.acfOriginal), ...
				'-', "DisplayName","Original ACF"); 
			plot(tc.lags, abs(acf), ...
				'--', "DisplayName", "Analystical ACF");

			legend("Location","northeast");
			xlabel("lags (-)");
			ylabel("Real Amplitude");
			grid on; 
		end
	end

end