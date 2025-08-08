classdef Simulation < handle
	properties (Access = private)
		N
		Q
		deltaD
	end

	methods
		%%
		function obj = Simulation(params)
			arguments
				params.N
				params.Q
			end

			obj.N = params.N;
			obj.Q = params.Q;

			obj.deltaD = obj.ComputeDeltaD();
		end

		%%
		function a = CoordinateDescent(obj, a0, options)
			arguments
				obj 
				a0 
				options.fftOversampling = 1
			end
			a = a0;
			
			for p = 3:obj.Q
				g = obj.ComputeG(a, p);

				G = fftshift(fft(g, options.fftOversampling*numel(g), 2));
				nG = numel(G);

				[~, k] = min(abs(G), [], "all");
				a(p) = 2*pi*(k-1)/nG;
			end
		end

		%%
		function isl_dB = ComputeIsl_dB(obj,a)
			paramPhase = zeros(1, obj.N);
			for n=1:obj.N
				paramPhase(n) = 0;
				for q=1:obj.Q
					paramPhase(n) = paramPhase(n) + ...
						a(q)*n^(q-1);
				end
			end

			x = exp(1i*sum(paramPhase, 1));
			
			acfFull = xcorr(x);
			acf = acfFull(obj.N:end);
			isl = 2*sum( abs(acf(2:end)).^2 );
			isl_dB = 10*log10(isl);
		end

	end

	methods (Access	= private)
		%%
		function deltaD = ComputeDeltaD(obj)
			deltaD = zeros(obj.N-1, obj.N-1, obj.N-1, obj.Q);
			for q = 1:obj.Q
				for k = 1:obj.N-1
					for n=1:obj.N-k
						for m=1:obj.N-k

							deltaD(n,m,k,q) = eqn.deltaD(n,m,k,q);

						end
					end
				end
			end
		end

		%%
		function g = ComputeG(obj, a, p)
			d_p = obj.deltaD(:,:,:,p);  
			dmin = min(d_p(:));  
			dmax = max(d_p(:));  
			d = dmin:dmax;
			id_offset = 1 - dmin;

			g = zeros(size(d));       
			for k = 1:obj.N-1
				for n = 1:obj.N-k
					for m = 1:obj.N-k

						id = obj.deltaD(n,m,k,p);
						idx = id + id_offset;

						phase = 0;
						for q = 2:obj.Q
							if q == p, continue; end
							phase = phase + a(q)*obj.deltaD(n,m,k,q);
						end
						g(idx) = g(idx) + 2*exp(1j*phase);

					end
				end
			end
		end

	end
end