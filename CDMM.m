classdef CDMM < handle
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
        save = @(varargin) [] % Initially NOP

        tmp_a_mm
        tmp_p
        tmp_iter_mm
	end

	methods
		function obj = CDMM(params)
			arguments
				params.nr_samples       (1,1) double {mustBePositive}           % Sequence Length
				params.nr_params        (1,1) double {mustBePositive}           % Number of Parameters
                params.nr_neighbors     (1,1) double {mustBePositive}   = 6    % Number of neighbors to include for interpolation
                params.nr_candidates    (1,1) double {mustBePositive}   = 30    % Number of minima to pursue
                params.max_mm_iters     (1,1) double {mustBePositive}   = 1e4   % Maximum number of MM steps
                params.nr_cd_iterations (1,1) double {mustBePositive}   = 1     % Number of complete CD iterations
                params.oversampling     (1,1) double {mustBePositive}   = 3     % Oversampling for NUFFT

                params.save_data        (1,1) logical                   = false
                params.file_location    (1,1) string {mustBeFolder}     = "02_Data"
                params.file_name        (1,1) string                    = NaN
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
            if params.save_data
                obj.data_logger = matfile(params.file_location + "\" + params.file_name + ".mat", "Writable",true);
                obj.save = @obj.update_log_file;

                % Prepare editable data
                obj.data_logger.params          = params;
                obj.data_logger.a               = zeros(obj.nr_params, 0);
                obj.data_logger.coordinate      = zeros(1,0);
                obj.data_logger.mm_step_count   = zeros(1,0);
                obj.data_logger.a_mm            = zeros(obj.nr_candidates, 0);
            end
		end

		function a = simulate(obj, a)
			arguments
				obj 
				a (:,1) double 
            end
			
            obj.save(a);

            for i_cd = 1:obj.nr_cd_iters
			    for p = 3:obj.nr_params
                    % Bring into DFT format
				    [g, d]      = eqn.G(a, p, obj.deltaD, obj.nr_samples, obj.nr_params); 
                    g           = g / ( obj.nr_samples^3 );
                    smoothness  = sum( d.^2 .* abs(g) ); % M-Smooth Lipschwitz continuous
                   
                    % Rescaling
                    energy  = sum( abs(g).^2 ); 
                    g       = g / sqrt(energy); 
    
                    % NUFFT for ISL and dISL/da
                    nufft_f = NUFFT("x",g, "K",obj.oversampling*numel(g));
                    nufft_f.prepareInterpolator("J",obj.nr_neighbors);
    
                    nufft_df = NUFFT("x",1i*d.*g, "K",obj.oversampling*numel(g));
                    nufft_df.prepareInterpolator("J",obj.nr_neighbors);
    
                    % Select multiple MM candidates
                    [~, idx_a] = findpeaks(- abs(nufft_f.X), "NPeaks", obj.nr_candidates-1, "SortStr","descend");
                    % [~, idx_a]  = mink( nufft_f.X , obj.nr_candidates-1 );
                    a_MM        = [ a(p) ; real(2*pi*nufft_f.f(idx_a)) ];
    
                    % Majorize / Minimize
                    for i = 1:obj.max_mm_iters
                        tangent = real(nufft_df.interpolate(a_MM));
                        % error = abs(nufft_df.error(a_MM)); % Interpolation error
                        % if abs(tangent) > error % Correction still save
                        % a_MM = a_MM - (tangent + sign(tangent).*error)/smoothness;
                        % end
                        a_MM = a_MM - tangent/smoothness;

                        obj.save(a_MM, p, i);
                    end
    
                    % Select final minima
                    isl     = sqrt(energy) * real(nufft_f.interpolate(a_MM));
				    [~, k]  = min(isl, [], "all");
                    a(p)    = a_MM(k);
                    obj.save(a);
                    clear a_MM phi
                end
            end
            obj.save();
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

    methods (Access=private)
        function update_log_file(obj, a, p, iter_mm)
            arguments
                obj 
                a       (:,1) double = NaN
                p       (1,1) double = NaN
                iter_mm (1,1) double = NaN
            end
            
            
            if ~isnan(a)
                if any([isnan(p), isnan(iter_mm)])
                    obj.data_logger.a(:,end+1) = a;
                    return;
                end
                obj.tmp_a_mm(:,end+1) = a;
                obj.tmp_p(end+1) = p;
                obj.tmp_iter_mm(end+1) = iter_mm;
            end

            len = size(obj.tmp_a_mm,2);
            if any([all([isnan(a); len > 0]); len > 99])
                obj.data_logger.a_mm(:,(end+1):(end+len)) = obj.tmp_a_mm;
                obj.data_logger.coordinate(:,(end+1):(end+len)) = obj.tmp_p;
                obj.data_logger.mm_step_count(:,(end+1):(end+len)) = obj.tmp_iter_mm;

                obj.tmp_a_mm = [];
                obj.tmp_p = [];
                obj.tmp_iter_mm = [];
                return;
            end
        end
    end
end