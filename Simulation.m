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
                
                g = g/(obj.N^3);
               
                E = sum(abs(g).^2); 
                g = g/sqrt(E); 

                nufft_f = NUFFT("x",g, "K",2*numel(g));
                nufft_df = NUFFT("x",1i*d.*g, "K",2*numel(g));

                J = 6;

                % Test Correctnes of the NUFFT for ISL
                % plot(nufft_f.f, 10*log10(real(nufft_f.X)*sqrt(E)),'-', "DisplayName","NUFFT");
                % for i = 1:numel(nufft_f.f)
                %     a0(p) = 2*pi*nufft_f.f(i);
                %     I(i) = obj.ComputeIsl_dB(a0);
                % end
                % hold on; 
                % plot(nufft_f.f, (I), '--', "DisplayName","ISL Calc");
                % legend("Location","best")
                

                nufft_f.prepareInterpolator("J",J);
                nufft_df.prepareInterpolator("J",J);

                % Find x-1 minimas
                [~, idx_a] = mink((nufft_f.X),29);

                a_MM(:,1) = [a(p); real(2*pi*nufft_f.f(idx_a))];
                phi(:,1) = real(nufft_f.interpolate(a_MM));

                M = sum(abs(d.^2.*g));
                Iter = 1e4;
                for iter = 1:Iter
                    a_MM(:,iter+1) = a_MM(:,iter) - (real(nufft_df.interpolate(a_MM(:,iter))) + real(nufft_df.error(a_MM(:,iter))))/(M);
                    phi(:,iter+1) = sqrt(E) * real(nufft_f.interpolate(a_MM(:,iter+1)));
                end
                % plot(a_MM(:,end)/(2*pi), 10*log10(phi(:,end)), 'o', "DisplayName","Optimized", "HandleVisibility","off");

				[isl, k] = min(phi(:,end), [], "all");
                if isl < obj.isl
                    obj.isl = isl;
                    a(p) = a_MM(k,end);
                    fprintf("+")
                else
                    fprintf("*");
                end
                clear a_MM phi
                
				
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