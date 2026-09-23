classdef General < handle

    properties (Access = private)
		params Algorithm.Params

		d
        g
        smoothness
        tol
        d_range
        periods
   	end

	methods
		function obj = General(params)
			arguments
				params Algorithm.Params
            end
			obj.params = params;

            [obj.d, obj.g] = calculate_gd(obj.params.nr_samples, obj.params.nr_params);

            obj.smoothness  = obj.calculate_smoothness();
            obj.d_range = obj.calculate_d_range();
            obj.periods = 2*pi./obj.calculate_gcd();
        end

        function a = run(obj, a, options)
			arguments
				obj 
				a (:,1) double 
                options.data_logger Datalogger.Type = Datalogger.Empty
            end
            change_dB_min = 0.5;
            delta_min = 10^(change_dB_min/20);
            f_tau = obj.params.tol + sqrt(obj.params.tol + 2*delta_min); 

            logger = options.data_logger;
            logger.log("params", obj.params);

            phase = a(3:end).'*obj.d;
            g_p = obj.g.*exp(1i*phase); 
			
            ISL_old = Inf;
            for i_cd = 1:obj.params.nr_cd_iterations
			    for p = 3:obj.params.nr_params
                    % Bring into DFT format
                    d_p = obj.d(p-2,:);
                    phase = a(p).'*d_p;
                    g_p = g_p.*exp(-1i*phase);

                    % Configure Search Space
                    tmp = numel(g_p) * abs(log10(1/obj.params.tol));
                    L = tmp / log10(tmp);

                    dA = pi/max(d_p);
                    if dA*L > obj.periods(p-2), L = obj.periods(p-2)/dA; end % Rescale L to one period
                    L = 2*ceil(L/2) - 1;  % Make odd; 

                    a_range = a(p) + dA*[-(L-1)/2;(L-1)/2];

                    % Configure Basins
                    if obj.params.nr_candidates == 1
                        a_grid = a(p); 
                    else
                        a_grid = linspace(a_range(1),a_range(2), obj.params.nr_candidates);
                        a_grid((obj.params.nr_candidates +1)/2) = a(p);
                    end

                    % NUFFT for ISL and dISL/da
                    f = oNufft.Type3(obj.d_range(p-2,:), a_range, "epsilon",obj.params.tol);
                    df = oNufft.Type3(obj.d_range(p-2,:), a_range, "epsilon",obj.params.tol);
                    f.apply(d_p, g_p);
                    df.apply(d_p, 1i*d_p.*g_p);
                    
                    % Select MM candidates
                    [~, idx_a] = min(real(f.interpolate(a_grid)));
                    a_mm = a_grid(idx_a);

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

                        if abs(tangent)  / sqrt(obj.smoothness(p-2)) < f_tau
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

                ISL_new = real(f.interpolate(a(p)));
                change_db = abs(20*log10(ISL_old) - 20*log10(ISL_new));
                if change_db < change_dB_min, break; end
                ISL_old = ISL_new; 
                
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

        function g = calculate_gcd(obj)
            g = obj.d(:,2); 
            for i = 4:2:size(obj.d,2)
                g = gcd(g, obj.d(:,i));
            end
        end
    end
end