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
                options.a_range = NaN
            end
            logger = options.data_logger;
            logger.log("params", obj.params);

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
                    if ~isnan(options.a_range)
                        a_grid = options.a_range(1):dA:options.a_range(end);
                        a_range = options.a_range;
                    else
                        tmp = obj.params.nr_samples^3 * log(obj.params.tol);
                        N = 2^ceil(log2(abs(tmp) / log(abs(tmp))));
                        %clear tmp;

                        f = unique(factor(p-1));
                        if isscalar(f)
                            period = 2*pi/f;
                        else
                            period = 2*pi;
                        end
                        
                        if N*dA > period
                            N = 2^ceil(log2(pi/dA));
                        end
                        a_grid = a(p) + dA*(-N/2:N/2).';
                        a_range = [a_grid(1), a_grid(end)];
                    end

                    f = oNufft.Type3(obj.d_range(p-2,:), a_range, "epsilon",obj.params.tol/40);
                    df = oNufft.Type3(obj.d_range(p-2,:), a_range, "epsilon",obj.params.tol/40);
                    f.apply(d_p, g_p);
                    df.apply(d_p, 1i*d_p.*g_p);
                    
                    % Select multiple MM candidates
                    a_mm = obj.select_mm_candidates(f, a_grid);

                    logger.log( ...
                        "isl",  f, ...
                        "a",    a, ...
                        "a_MM", a_mm, ...
                        "i_MM", 0, ...
                        "p",    p, ...
                        "i_cd", i_cd ...
                        );

                    % Majorize / Minimize
                    for i = 1:obj.params.max_mm_iters
                        tangent = real(df.interpolate(a_mm));
                        a_step = - tangent/obj.smoothness(p-2);
                        
                        if ~allbetween(a_mm + a_step, a_range(1), a_range(2))
                            break; % Out of reach for save NUFFT interpolation
                        end

                        if all([abs(tangent / obj.smoothness(p-2)) < obj.params.tol, i > 1])
                            break; % No more progress
                        end

                        a_mm = a_mm + a_step;
                        logger.log( ...
                            "a_MM", a_mm, ...
                            "i_MM", i, ...
                            "p",    p, ...
                            "i_cd", i_cd ...
                            );
                    end
    
                    % Select final minima
                    isl     = real(f.interpolate(a_mm));
				    [~, k]  = min(isl);
                    a(p)    = a_mm(k);

                    phase = a(p).'*d_p;
                    g_p = g_p.*exp(1i*phase);
                end
            end

            logger.log("a", a);
            %delete(logger);
        end

    end

    methods (Access = private)
        function [a_mm] = select_mm_candidates(obj, f, candidates)
            arguments
                obj 
                f oNufft.Type3
                candidates 
            end

            if obj.params.nr_candidates > 1
                [~, idx_a] = findpeaks(- real(f.interpolate(candidates)), "NPeaks", obj.params.nr_candidates, "SortStr","descend");
                if numel(idx_a) < obj.params.nr_candidates
                    warning("CDMM: Could not find the required nr of MM candidates");
                end
                a_mm        = [ candidates(idx_a) ];
            else
                [~, idx_a] = min(real(f.interpolate(candidates)));
                a_mm = candidates(idx_a);
            end
        end

        function smoothness = calculate_smoothness(obj)
            [~, I] = sort(obj.d(1,:)); % Sorting can be done in any parameter dimension
            d_ = obj.d(:,I);
            g_ = obj.g(I);

            N = floor(numel(g_)/2);
            g_(1:N) = [];
            d_(:,1:N) = [];
            smoothness  = sum( d_.^2 .* g_ ,2); % M-Smooth Lipschwitz continuous
        end

        function d_ranges = calculate_d_range(obj)
            d_ranges = [min(obj.d, [], 2), max(obj.d, [], 2)];
        end
    end
end