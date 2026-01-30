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
        log_descent    = @(varargin) [] % Initially NOP
        log_mm_step    = @(varargin) [] % Initially NOP
        save_all        = @(varargin) [] % Initially NOP
            
        tmp_a_mm
        tmp_nufft
        tmp_a
	end

	methods
		function obj = CDMM(params)
			arguments
				params.nr_samples       (1,1) double {mustBePositive}           % Sequence Length
				params.nr_params        (1,1) double {mustBePositive}           % Number of Parameters
                params.nr_neighbors     (1,1) double {mustBePositive}   = 4    % Number of neighbors to include for interpolation
                params.nr_candidates    (1,1) double {mustBePositive}   = 30    % Number of minima to pursue
                params.max_mm_iters     (1,1) double {mustBePositive}   = 1e4   % Maximum number of MM steps
                params.nr_cd_iterations (1,1) double {mustBePositive}   = 1     % Number of complete CD iterations
                params.oversampling     (1,1) double {mustBePositive}   = 5     % Oversampling for NUFFT

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
                obj.log_descent = @obj.update_descent;
                obj.log_mm_step = @obj.update_a_mm;
                obj.save_all    = @obj.save_file;

                % Prepare editable data
                obj.data_logger.params = params;
            end
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
                   
                    % Rescaling
                    % energy  = sum( abs(g).^2 ); 
                    % g       = g / sqrt(energy); 
    
                    % NUFFT for ISL and dISL/da
                    nufft_f = NUFFT("x",g, "K",obj.oversampling*numel(g));
                    nufft_f.prepareInterpolator("J",obj.nr_neighbors);
    
                    nufft_df = NUFFT("x",1i*d.*g, "K",obj.oversampling*numel(g));
                    nufft_df.prepareInterpolator("J",obj.nr_neighbors);

                    obj.log_descent(nufft_f, a, p, i_cd);
    
                    % Select multiple MM candidates
                    [~, idx_a] = findpeaks(- abs(nufft_f.X), "NPeaks", obj.nr_candidates-1, "SortStr","descend");
                    % [~, idx_a]  = mink( nufft_f.X , obj.nr_candidates-1 );
                    a_MM        = [ a(p) ; real(2*pi*nufft_f.f(idx_a)) ];
    
                    obj.log_mm_step(a_MM, 0, p, i_cd);
                    % Majorize / Minimize
                    for i = 1:obj.max_mm_iters
                        tangent = real(nufft_df.interpolate(a_MM));
                        % error = abs(nufft_df.error(a_MM)); % Interpolation error
                        % if abs(tangent) > error % Correction still save
                        % a_MM = a_MM - (tangent + sign(tangent).*error)/smoothness;
                        % end
                        a_MM = a_MM - tangent/smoothness;
                        a_MM = mod(a_MM, 2*pi);
                        obj.log_mm_step(a_MM, i, p, i_cd);
                    end
    
                    % Select final minima
                    isl     = real(nufft_f.interpolate(a_MM));
				    [~, k]  = min(isl, [], "all");
                    a(p)    = a_MM(k);
                    clear a_MM phi
                end
            end

            obj.log_descent(nufft_f, a, p, i_cd);
            obj.save_all();
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
        function update_a_mm(obj, a, i, p, i_cd)
            arguments
                obj 
                a       (:,1) double = NaN
                i       (1,1) double = NaN
                p       (1,1) double = NaN
                i_cd    (1,1) double = NaN
            end
            
            obj.tmp_a_mm(:,i+1,p-2,i_cd) = a;
        end

        function update_descent(obj, nufft, a, p, i_cd)
            obj.tmp_nufft{p-2, i_cd} = nufft;
            obj.tmp_a(:, p-2, i_cd) = a;
        end

        function save_file(obj)
            obj.data_logger.a_mm = obj.tmp_a_mm;
            obj.data_logger.a = obj.tmp_a;
            obj.data_logger.nufft = obj.tmp_nufft;
        end
    end
end