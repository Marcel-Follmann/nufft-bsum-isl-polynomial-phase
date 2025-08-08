classdef TestMathematicalDerivations < matlab.unittest.TestCase
	properties
		N = 3
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
				for q = 1:tc.Q
					for k = 1:tc.N-1
						for n=1:tc.N-k
							for m=1:tc.N-k
								deltaD(n,m,k,q) = eqn.deltaD(n,m,k,q);
							end
						end
					end
				end
				
				% Extract the p-th differences over all (n,m,k):
				dp = squeeze(deltaD(:,:,:,p));  
				% Compute min and max possible Δd^(p)
				dmin = min(dp(:));  
				dmax = max(dp(:));  
				dVals = dmin:dmax;            % vector of all possible integer bins
				offset = 1 - dmin;            % so that d = 0 maps to index = offset
				G = zeros(size(dVals));       % preallocate G_d

				for k = 1:(tc.N-1)
  
  					% — Handle the diagonal m=n contribution (Δd^(p)=0, A=1):
  					G(offset) = G(offset) + 2*(tc.N-k);
					
  					% — Now off-diagonal n<m only:
  					for n = 1:(tc.N-k)
    					for m = (n+1):(tc.N-k)
      					%  p-th difference:
      					d = deltaD(n,m,k,p);
      					idx = d + offset;        % bin index
					
      					%  phase sum over all q≠p:
      					phi = 0;
      					for q = 2:tc.Q
        					if q==p, continue; end
        					phi = phi + tc.a(q)*deltaD(n,m,k,q);
      					end
					
      					% accumulate 4⋅cos(φ):
      					G(idx) = G(idx) + 4*cos(phi);
    					end
					end

				end

				% --- after your loops that fill G and you still have dmin, dmax, offset ---

				% 1) set up zero-padded FFT length for fine angular resolution
				L       = length(G);                      % original # of bins
				M       = 4 * 2^nextpow2(L);              % e.g. 4× next power-of-2
				angles  = linspace(0, 2*pi, M);           % candidate a_p in [0,2π)
				
				% 2) embed G into a length-M vector H (so that d=0 → index offset)
				H = zeros(1, M);
				H(offset + (0:(L-1))) = G;
				
				% 3) FFT and find coarsest minimizer
				P       = fft(H);
				[~, kmin] = min(abs(P));                  % index of smallest magnitude
				a_p0    = angles(kmin);                   % grid-based estimate
				
				% 4) Optional: refine via fminbnd on the continuous cost
				dVals = dmin:dmax;  
				cost  = @(ap) abs( sum( G .* exp(1j * ap * dVals) ) );
				delta = 2*pi/M;                            % half-spacing of your grid
				a_lo  = a_p0 - delta;
				a_hi  = a_p0 + delta;
				opts  = optimset('TolX',1e-10,'Display','off');
				[a_p_opt, fval] = fminbnd(cost, a_lo, a_hi, opts);
				
				% 5) display results
				fprintf('FFT‐grid estimate:   a_p = %.8f rad\n', a_p0);
				fprintf('Refined (fminbnd):   a_p = %.12f rad\n', a_p_opt);


				
			end
			Subtest_Eqn_8(tc,3);
		end

	end

end