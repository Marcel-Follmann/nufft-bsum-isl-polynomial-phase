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
                
                % g(d>0) = [];
                % d(d>0) = [];
                g = g/(obj.N^3);
               
                E = sum(g.^2); 
                g = g/sqrt(E);

                nufft_f = NUFFT("x",g, "K",4*numel(g));
                nufft_f2 = NUFFT("x",g, "K",10*numel(g));
                nufft_df = NUFFT("x",1i*d.*g, "K",4*numel(g));
                J = 10;
                nufft_f.prepareInterpolator("J",J);
                nufft_df.prepareInterpolator("J",J);

                % Find 5 minimas
                [~, idx_a] = mink((nufft_f.X),20);
                a_MM = real(nufft_f.f(idx_a));
                phi = real(nufft_f.interpolate(a_MM*2*pi));
                M = sum(d.^2.*g);
                Iter = 1e2;
                plot(nufft_f.f,10*log10(nufft_f.X*sqrt(E)), 'o', "DisplayName", 'Original FFT');hold on; 
                % plot(nufft_f2.f,10*log10(nufft_f2.X*sqrt(E)), '--', "DisplayName", 'Zero padded FFT');
                for iter = 1:Iter
                    a_MM(:,iter+1) = a_MM(:,iter) + (real(nufft_df.interpolate(a_MM(:,iter)*2*pi)) + real(nufft_df.error(a_MM(:,iter)*2*pi)))/(M);
                    phi(:,iter+1) = real(nufft_f.interpolate(a_MM(:,iter+1)*2*pi));
                    plot(a_MM(:,iter+1), 10*log10(phi(:,iter+1)*sqrt(E)), '.', "Color", iter/Iter*[1,0,0]); 
                end
                title("Red is newer")
                
                phiFft = fft(ifftshift(g), numel(g));
				nFft = numel(phiFft);

                close all; 

                nNufft = nFft;
                aNufft = 0:1/nFft:1/2-1/nNufft;

                phiNufft = nufft_f.interpolate(aNufft*(2*pi));


                nexttile; hold on; 
                plot(d, 2*abs(g), "DisplayName","FFT"); 
                legend("Location","northeast");

                nexttile; hold on; 
                plot(0:1/nFft:1-1/nFft, (phiFft), "DisplayName", "FFT");
                plot(aNufft, (phiNufft), ".", "DisplayName","NUFFT"); 
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