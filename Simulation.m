classdef Simulation < handle
	properties (Access = private)
		N
		Q
		deltaD
        isl
	end

	methods
		%%
		function obj = Simulation(params)
			arguments
				params.N (1,1) double {mustBePositive} % Sequence Length
				params.Q (1,1) double {mustBePositive} % Number of Parameters
			end

			obj.N = params.N;
			obj.Q = params.Q;
            obj.isl = Inf;

			obj.deltaD = eqn.deltaD_mat(obj.N, obj.Q);
            fprintf("Initilization complete.\n")
		end

		%%
		function a = CoordinateDescent(obj, a0, options)
			arguments
				obj 
				a0 (:,1) double % Initial Parameters
				options.fftOversampling (1,1) double {mustBePositive} = 1
            end
            fprintf("Start CD step. ");
			a = a0;
			
			for p = 3:obj.Q
				g = eqn.G(a, p, obj.deltaD, obj.N, obj.Q);
                gShift = ifftshift(g);
				G = ifft(gShift, options.fftOversampling*numel(g));
				nG = numel(G);

				[isl, k] = min(abs(G), [], "all");
                if isl < obj.isl
                    obj.isl = isl;
                    a(p) = obj.fftBin2parameterVaule(k,nG);
                    fprintf("+")
                else
                    fprintf("*");
                end
                
				
            end
            fprintf(" Done.\n");
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
			isl = 2*sum( abs(acf(2:end)).^2 ./ obj.N^3);
			isl_dB = 10*log10(isl);
		end

	end

	methods (Static)
        %%
        function a = fftBin2parameterVaule(bin, nBins)
            a_unwrapped = 2*pi*(bin-1)/nBins;
            a = mod(a_unwrapped + pi, 2*pi) - pi;
        end

	end
end