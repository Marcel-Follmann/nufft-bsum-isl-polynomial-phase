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
		function a = CoordinateDescent(obj, a0)
			arguments
				obj 
				a0 (:,1) double % Initial Parameters
            end
            fprintf("Start CD step. ");
			a = a0;
			
			for p = 3:obj.Q
				[g, d] = eqn.G(a, p, obj.deltaD, obj.N, obj.Q);
                g = g/(obj.N^2);
                phiFft = ifft(ifftshift(g), numel(g));
				nFft = numel(phiFft);
                % g(g<128) = 0;
                dNufft = d;
                gNufft = g; 

                dNufft(abs(g) == 0) = [];
                gNufft(abs(g) == 0) = [];

                nNufft = .010*nFft;
                aNufft = 0:1/nFft:1/2-1/nNufft;
                % atest = mod(atest + pi, 2*pi) - pi;

                phiNufft = nufft(flip(gNufft), dNufft, aNufft) / nFft;

                nexttile; hold on; 
                plot(d, 2*abs(g), "DisplayName","FFT"); 
                plot(dNufft, 2*abs(gNufft), ".", "DisplayName","NUFFT");
                legend("Location","northeast");

                nexttile; hold on; 
                plot(0:1/nFft:1-1/nFft, abs(phiFft), "DisplayName", "FFT");
                plot(aNufft, abs(phiNufft), ".", "DisplayName","NUFFT"); 
                legend("Location","northeast")

                fprintf("Energy: %.2e\n", sum(abs(g).^2));

				[isl, k] = min(abs(phiFft), [], "all");
                if isl < obj.isl
                    obj.isl = isl;
                    a(p) = obj.fftBin2parameterVaule(k,nFft);
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