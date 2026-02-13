classdef CDMMT3 < handle
	properties (Access = private)
		nr_samples
		nr_params
        nr_neighbors
        nr_candidates
        nr_cd_iters
        max_mm_iters
        oversampling
		dg
        Gd
        smoothness
        tol

        a_minima

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
                params.tol              (1,1) double {mustBeNonnegative}   = 0  % Tolerance for local minima

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
            obj.tol             = params.tol;

            % Prepare DFT reformatting
			% obj.deltaD = eqn.deltaD_mat(obj.nr_samples, obj.nr_params);

            [obj.dg, Gd_abs] = calculate_gd(obj.nr_samples, obj.nr_params);
            obj.Gd = 2*Gd_abs./(obj.nr_samples^3);
            [~, I] = sort(obj.dg(1,:));
            d = obj.dg(:,I);
            g = obj.Gd(I);


            N = floor(numel(g)/2);
            g(1:N) = [];
            d(:,1:N) = [];
            d_max = max(d, [], 'all');
            for i = 1:size(d,1)
                tmp(i,:) = accumarray(d(i,:).'+1, g.', [d_max+1,1]);
            end
            
            n_tmp = 0:size(tmp,2)-1;
            % stem(tmp(2:end)); hold on;
            obj.smoothness  = sum( n_tmp.^2 .* tmp ,2); % M-Smooth Lipschwitz continuous
            
            a_minima_mean = 2*pi./sqrt(sum(n_tmp.^2 .* tmp.^2, 2)./sum(tmp.^2,2));
            p = 0.99; % 99th percentile
            obj.a_minima = -a_minima_mean*log(1-p);


            % Prepare data logger
            obj.data_logger = DataLogger("log_enabled",params.save_data, "file_location",params.file_location, "file_name",params.file_name);
            obj.data_logger.log("params", params);
		end

		function a = type2(obj, a)
			arguments
				obj 
				a (:,1) double 
            end

            phase = a(3:end).'*obj.dg;
            g = obj.Gd.*exp(-1i*phase); % 2*(2, because of Real part)
			
            for i_cd = 1:obj.nr_cd_iters
			    for p = 3:obj.nr_params
                    % Bring into DFT format
                    d = obj.dg(p-2,:);
                    phase = a(p).'*d;
                    g = g.*exp(1i*phase);

                    dA = 2*pi/max(d);

                    % NUFFT for ISL and dISL/da
                    % candidates(:,1) = linspace(0, 2*pi, 1e4);
                    % candidates = a(p) + dA*(-obj.nr_candidates*n:obj.nr_candidates*n).';
                    candidates = a(p) + dA*(-obj.nr_samples^2:obj.nr_samples^2).';
                    a_range = [candidates(1), candidates(end)];
                    f = oNufft.Type3(d, candidates, "epsilon",1e-16);
                    df = oNufft.Type3(d, candidates, "epsilon",1e-16);
                    f.apply(d, g);
                    df.apply(d, 1i*d.*g);
                    
                    obj.data_logger.log("isl", f, "a", a);
    
                    % Select multiple MM candidates
                    if obj.nr_candidates > 1
                        [~, idx_a] = findpeaks(- real(f.interpolate(candidates)), "NPeaks", obj.nr_candidates, "SortStr","descend");
                        if numel(idx_a) < obj.nr_candidates
                            warning("CDMM: Could not find the required nr of MM candidates");
                        end
                        a_MM        = [ candidates(idx_a) ];
                    else
                        a_MM = a(p);
                    end
    
                    obj.data_logger.log("a_MM", a_MM, "i_MM", 0, "p", p, "i_cd", i_cd);
                    % Majorize / Minimize
                    for i = 1:obj.max_mm_iters
                        tangent = -real(df.interpolate(a_MM));
                        
                        a_MM = a_MM - tangent/obj.smoothness(p-2);
                        if allbetween(a_MM, a_range(1), a_range(2))
                            obj.data_logger.log("a_MM", a_MM, "i_MM", i, "p", p, "i_cd", i_cd);
                            if all(abs(tangent / obj.smoothness(p-2)) < obj.tol)
                                break;
                            end
                        else
                            break;
                        end
                    end
    
                    % Select final minima
                    isl     = real(f.interpolate(a_MM));
				    [~, k]  = min(isl);
                    a(p)    = a_MM(k);
                    % clear a_MM
                    phase = a(p).'*d;
                    g = g.*exp(-1i*phase);
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