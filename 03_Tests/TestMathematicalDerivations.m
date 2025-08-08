classdef TestMathematicalDerivations < matlab.unittest.TestCase
	properties
		N = 20
		Q = 4
		a
		x
		lags
		acfOriginal
		islOriginal
		relTol = 1e-3
	end

	methods (TestClassSetup)
		function SetupClass(tc)
			close all; 
			
			tc.a = randn(tc.Q, 1);
			
			paramPhase = zeros(1, tc.N);
			for n=1:tc.N
				paramPhase(n) = 0;
				for q=1:tc.Q
					paramPhase(n) = paramPhase(n) + ...
						tc.a(q)*n^(q-1);
				end
			end

			tc.x = exp(1i*sum(paramPhase, 1));
			
			tc.lags = 0:(tc.N-1);
			
			acfOriginal_full = xcorr(tc.x);
			tc.acfOriginal = acfOriginal_full(tc.N:end);
			tc.islOriginal = 2*sum( abs(tc.acfOriginal(2:end)).^2 );
		end
	end

	methods (TestMethodSetup)
		
	end
		
	methods (Test)
		% Test methods

		function Test_Eqn_2(tc) 
			% r_xx(k) = \sum_{n=1}^{N-k} x_n^* * x_{n+k}, k = {0, 1, ..., N-1}
			
			acfAnalytical = zeros(size(tc.acfOriginal)); 
			for k = 0:tc.N-1
			
				for n = 1:tc.N-k
					acfAnalytical(k+1) = acfAnalytical(k+1) + ...
						conj(tc.x(n)) * tc.x(n+k);
				end
			end


			tc.verifyEqual(abs(acfAnalytical), abs(tc.acfOriginal), "RelTol", tc.relTol);
		end

		function Test_Eqn_4(tc) 
			acfAnalytical = zeros(size(tc.acfOriginal));
			for k = 0:tc.N-1
				for n = 1:tc.N-k
					phaseSum = 0;
					for q=2:tc.Q 
						phaseSum = phaseSum + ...
							tc.a(q) * eqn.d(n,k,q);
					end

					acfAnalytical(k+1) = acfAnalytical(k+1) + ...
						exp(1j*phaseSum);
				end
			end

			tc.verifyEqual(acfAnalytical, tc.acfOriginal, "RelTol", tc.relTol);
		end

		function Test_Eqn_5(tc) 
			function Subtest_Eqn_5(tc,p)
				acfAnalytical = zeros(size(tc.acfOriginal));
				for k = 0:tc.N-1
					
					for n = 1:tc.N-k
						acfAnalytical(k+1) = acfAnalytical(k+1) + ...
							exp( 1j* tc.a(p) * eqn.d(n,k,p) ) * eqn.c(n,k,tc.a, p);
					end
				end
			
				tc.verifyEqual(abs(acfAnalytical), abs(tc.acfOriginal), ...
					sprintf("Failed for p = %d", p), ...
					"RelTol", tc.relTol);
			end

			for p = 1:tc.Q
				Subtest_Eqn_5(tc,p);
			end
		end

		function Test_Eqn_7(tc)
			function Subtest_Eqn_7(tc, p)
				islAnalytical = 0;
				for k = 1:tc.N-1
					
					for n=1:tc.N-k
						for m=1:tc.N-k
							islAnalytical = islAnalytical + ...
								2 * eqn.A(n,m,k,tc.a,p) * ...
								exp( 1j*tc.a(p) * eqn.deltaD(n,m,k,p) );
						end
					end
					
				end
	
				tc.verifyEqual(islAnalytical, tc.islOriginal, ...
					sprintf("Failed for p = %d", p), ...
					"RelTol", tc.relTol);
			end

			for p = 1:tc.Q
				Subtest_Eqn_7(tc,p);
			end
		end

		function Test_Eqn_8(tc)
			function Subtest_Eqn_8(tc,p)
				% Compoute all possible deltaD
				for q = 1:tc.Q
					for k = 1:tc.N-1
						for n=1:tc.N-k
							for m=1:tc.N-k
								deltaD(n,m,k,q) = eqn.deltaD(n,m,k,q);
							end
						end
					end
				end
				
				% Construct Gd
				d_p = deltaD(:,:,:,p);  
				dmin = min(d_p(:));  
				dmax = max(d_p(:));  
				d = dmin:dmax;
				id_offset = 1 - dmin;

				G = zeros(size(d));       
				for k = 1:(tc.N-1)
					for n = 1:tc.N-k
						for m = 1:tc.N-k
							id = deltaD(n,m,k,p);
							idx = id+id_offset;

							phase = 0;
							for q = 2:tc.Q
								if q == p, continue; end
								phase = phase + tc.a(q)*deltaD(n,m,k,q);
							end
							G(idx) = G(idx) + 2*exp(1j*phase);
						end
					end
				end

				islAnalytical = 0;
				for id = 1:numel(d)
					islAnalytical = islAnalytical + ...
						G(id)*exp(1j*tc.a(p)*d(id));
				end

				tc.verifyEqual(islAnalytical, tc.islOriginal, ...
					sprintf("Failed for p = %d", p), ...
					"RelTol", tc.relTol);
				
			end
			
			for p = 1:tc.Q
				Subtest_Eqn_8(tc,p);
			end
		end

	end

end