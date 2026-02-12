classdef CDMMT3 < handle
	properties (Access = private)
		nr_samples
		nr_params
        nr_neighbors
        nr_candidates
        nr_cd_iters
        max_mm_iters
        oversampling
		deltaD

        data_logger
   	end

	methods
		function obj = CDMMT3(params)
			arguments
				params.nr_samples       (1,1) double {mustBePositive}           % Sequence Length
				params.nr_params        (1,1) double {mustBePositive}           % Number of Parameters
                params.nr_neighbors     (1,1) double {mustBePositive}   = 4    % Number of neighbors to include for interpolation
                params.nr_candidates    (1,1) double {mustBePositive}   = 30    % Number of minima to pursue
                params.max_mm_iters     (1,1) double {mustBePositive}   = 1e3   % Maximum number of MM steps
                params.nr_cd_iterations (1,1) double {mustBePositive}   = 1     % Number of complete CD iterations
                params.oversampling     (1,1) double {mustBePositive}   = 2     % Oversampling for NUFFT

                params.save_data        (1,1) logical                   = false
                params.file_location    (1,1) string {mustBeFolder}     = "02_Data"
                params.file_name        (1,1) string                    = ""
			end

			obj.nr_samples      = params.nr_samples;
			obj.nr_params       = params.nr_params;
            obj.nr_neighbors    = params.nr_neighbors;
            obj.max_mm_iters    = params.max_mm_iters;
            obj.nr_candidates   = params.nr_candidates;
            obj.nr_cd_iters     = params.nr_cd_iterations;
            obj.oversampling    = params.oversampling;

            % Prepare DFT reformatting
			obj.deltaD = eqn.deltaD_mat(obj.nr_samples, obj.nr_params);

            % Prepare data logger
            obj.data_logger = DataLogger("log_enabled",params.save_data, "file_location",params.file_location, "file_name",params.file_name);
            obj.data_logger.log("params", params);
		end

		function a = type2(obj, a)
			arguments
				obj 
				a (:,1) double 
            end
			
            for i_cd = 1:obj.nr_cd_iters
			    for p = 3:obj.nr_params
                    % Bring into DFT format
				    [g, d]      = eqn.G(a, p, obj.deltaD, obj.nr_samples, obj.nr_params); 
                    g           = g / ( obj.nr_samples^3 );
                    smoothness  = sum( d.^2 .* abs(g) ); % M-Smooth Lipschwitz continuous
                   
                    d(g == 0) = [];
                    g(g == 0) = [];
                    % g(d<0) = [];
                    % d(d<0) = [];

                    dA = 20* 2*pi/max(d);

                    % NUFFT for ISL and dISL/da
                    candidates(:,1) = a(p) + 0.5*dA*linspace(-1, 1, obj.nr_candidates*100+1);
                    f = oNufft.Type3(d, candidates);
                    df = oNufft.Type3(d, candidates);
                    f.apply(d, g);
                    df.apply(d, 1i*d.*g);

                    obj.data_logger.log("isl", f, "a", a);
    
                    % Select multiple MM candidates
                    if obj.nr_candidates > 1
                        [~, idx_a] = findpeaks(- abs(f.interpolate(candidates)), "NPeaks", obj.nr_candidates, "SortStr","descend");
                        a_MM        = [ candidates(idx_a) ];
                    else
                        a_MM = a(p);
                    end
    
                    obj.data_logger.log("a_MM", a_MM, "i_MM", 0, "p", p, "i_cd", i_cd);
                    % Majorize / Minimize
                    for i = 1:obj.max_mm_iters
                        tangent = -real(df.interpolate(a_MM));
                        if all(tangent / smoothness < 1e-5)
                            break;
                        end
                        a_MM = a_MM - tangent/smoothness;
                        % a_MM = mod(a_MM, 2*pi);
                        obj.data_logger.log("a_MM", a_MM, "i_MM", i, "p", p, "i_cd", i_cd);
                    end
    
                    % Select final minima
                    isl     = real(f.interpolate(a_MM));
				    [~, k]  = min(isl, [], "all");
                    a(p)    = a_MM(k);
                    % clear a_MM phi
                end
            end

            obj.data_logger.log("a", a);
        end

        function delete(obj)
            delete(obj.data_logger);
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