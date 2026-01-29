classdef CdNufftMm < handle
	properties (Access = private)
		nr_samples
		nr_params
        nr_neighbors
        nr_candidates
        max_iters
		deltaD
	end

	methods
		function obj = CdNufftMm(params)
			arguments
				params.nr_samples       (1,1) double {mustBePositive}       % Sequence Length
				params.nr_params        (1,1) double {mustBePositive}       % Number of Parameters
                params.nr_neighbors     (1,1) double {mustBePositive} = 6   % Number of neighbors to include for interpolation
                params.nr_candidates    (1,1) double {mustBePositive} = 30  % Number of minima to pursue
                params.max_iters        (1,1) double {mustBePositive} = 1e4 % Maximum number of MM steps
			end

			obj.nr_samples      = params.nr_samples;
			obj.nr_params       = params.nr_params;
            obj.nr_neighbors    = params.nr_neighbors;
            obj.max_iters       = params.max_iters;
            obj.nr_candidates   = params.nr_candidates;

            % Prepare DFT reformatting
			obj.deltaD = eqn.deltaD_mat(obj.nr_samples, obj.nr_params);
		end

		function a = CoordinateDescent(obj, a0)
			arguments
				obj 
				a0 (:,1) double 
            end
			a = a0; 
			
			for p = 3:obj.nr_params
                % Bring into DFT format
				[g, d]      = eqn.G(a, p, obj.deltaD, obj.nr_samples, obj.nr_params); 
                g           = g / ( obj.nr_samples^3 );
                smoothness  = sum( d.^2 .* abs(g) ); % M-Smooth Lipschwitz continuous
               
                % Rescaling
                energy  = sum( abs(g).^2 ); 
                g       = g / sqrt(energy); 

                % NUFFT for ISL and dISL/da
                nufft_f = NUFFT("x",g, "K",2*numel(g));
                nufft_f.prepareInterpolator("J",obj.nr_neighbors);

                nufft_df = NUFFT("x",1i*d.*g, "K",2*numel(g));
                nufft_df.prepareInterpolator("J",obj.nr_neighbors);

                % Select multiple MM candidates
                [~, idx_a]  = mink( nufft_f.X , obj.nr_candidates-1 );
                a_MM        = [ a(p) ; real(2*pi*nufft_f.f(idx_a)) ];

                % Majorize / Minimize
                for i = 1:obj.max_iters
                    tangent = real(nufft_df.interpolate(a_MM)) + ...
                        real(nufft_df.error(a_MM)); % Corrected for interpolation error
                    a_MM    = a_MM - tangent/smoothness;
                end

                % Select final minima
                isl     = sqrt(energy) * real(nufft_f.interpolate(a_MM));
				[~, k]  = min(isl, [], "all");
                a(p)    = a_MM(k);
                clear a_MM phi
            end
		end

		function isl_dB = ComputeIsl_dB(obj,a)
			paramPhase = zeros(1, obj.nr_samples);
			for n=1:obj.nr_samples
				paramPhase(n) = 0;
				for q=1:obj.nr_params
					paramPhase(n) = paramPhase(n) + ...
						a(q)*n^(q-1);
				end
			end

			x = exp(1i*sum(paramPhase, 1));
			
			acfFull = xcorr(x)/obj.nr_samples;
			acf = acfFull(obj.nr_samples:end);
			isl_lin = 2*sum( abs(acf(2:end)).^2)/obj.nr_samples;
			isl_dB = 10*log10(isl_lin);
		end

    end
end