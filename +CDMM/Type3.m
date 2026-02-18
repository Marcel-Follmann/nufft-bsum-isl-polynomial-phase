classdef Type3 < handle
	properties (Access = private)
		params Params

		d
        g
        smoothness
        tol
        d_range
   	end

	methods
		function obj = Type3(params)
			arguments
				params Params
            end
			obj.params = params;

            [obj.d, obj.g] = calculate_gd(obj.params.nr_samples, obj.params.nr_params);

            obj.smoothness  = obj.calculate_smoothness();
            obj.d_range = obj.calculate_d_range();
		end

        function a = run(obj, a, options)
			arguments
				obj 
				a (:,1) double 
                options.data_logger Datalogger.Type = Datalogger.Empty
            end
            options.data_logger.log("params", obj.params);

            phase = a(3:end).'*obj.d;
            g_p = obj.g.*exp(1i*phase); 
			
            for i_cd = 1:obj.params.nr_cd_iterations
			    for p = 3:obj.params.nr_params
                    % Bring into DFT format
                    d_p = obj.d(p-2,:);
                    phase = a(p).'*d_p;
                    g_p = g_p.*exp(-1i*phase);

                    

                    % NUFFT for ISL and dISL/da
                    dA = pi/max(d_p);
                    candidates = a(p) + dA*(-obj.params.nr_samples:obj.params.nr_samples).';
                    a_range = [candidates(1), candidates(end)];

                    f = oNufft.Type3(obj.d_range(p-2,:), a_range, "epsilon",obj.params.tol/1000);
                    df = oNufft.Type3(obj.d_range(p-2,:), a_range, "epsilon",obj.params.tol/1000);
                    f.apply(d_p, g_p);
                    df.apply(d_p, 1i*d_p.*g_p);
                    
                    options.data_logger.log("isl", f, "a", a);
    
                    % Select multiple MM candidates
                    if obj.params.nr_candidates > 1
                        [~, idx_a] = findpeaks(- real(f.interpolate(candidates)), "NPeaks", obj.params.nr_candidates, "SortStr","descend");
                        if numel(idx_a) < obj.params.nr_candidates
                            warning("CDMM: Could not find the required nr of MM candidates");
                        end
                        a_MM        = [ candidates(idx_a) ];
                    else
                        [~, idx_a] = min(real(f.interpolate(candidates)));
                        a_MM = candidates(idx_a);
                    end
    
                    options.data_logger.log("a_MM", a_MM, "i_MM", 0, "p", p, "i_cd", i_cd);

                    % Majorize / Minimize
                    for i = 1:obj.params.max_mm_iters
                        tangent = real(df.interpolate(a_MM));
                        
                        a_MM = a_MM - tangent/obj.smoothness(p-2);
                        if allbetween(a_MM, a_range(1), a_range(2))
                            options.data_logger.log("a_MM", a_MM, "i_MM", i, "p", p, "i_cd", i_cd);
                            if all(abs(tangent / obj.smoothness(p-2)) < obj.params.tol)
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

                    phase = a(p).'*d_p;
                    g_p = g_p.*exp(1i*phase);
                end
            end

            options.data_logger.log("a", a);
            delete(options.data_logger);
        end

        function smoothness = calculate_smoothness(obj)
            [~, I] = sort(obj.d(1,:)); % Sorting can be done in any parameter dimension
            d_ = obj.d(:,I);
            g_ = obj.g(I);

            N = floor(numel(g_)/2);
            g_(1:N) = [];
            d_(:,1:N) = [];
            d_max = max(d_, [], 'all');
            for i = 1:size(d_,1)
                g_accum(i,:) = accumarray(d_(i,:).'+1, g_.', [d_max+1,1]);
            end
            
            d_accum = 0:size(g_accum,2)-1;
            smoothness  = sum( d_accum.^2 .* g_accum ,2); % M-Smooth Lipschwitz continuous
        end

        function d_ranges = calculate_d_range(obj)
            d_ranges = [min(obj.d, [], 2), max(obj.d, [], 2)];
        end

    end
end