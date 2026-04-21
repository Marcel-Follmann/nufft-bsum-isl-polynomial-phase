classdef Quadratic < handle

    properties (Access = private)
		params Algorithm.Params

		d
        g
        smoothness
        tol
        d_range

        f, df
   	end

	methods
		function obj = Quadratic(params)
			arguments
				params Algorithm.Params
            end
			obj.params = params;

            [obj.d, obj.g] = calculate_gd(obj.params.nr_samples, 3);

            obj.smoothness  = obj.calculate_smoothness();
            obj.d_range = obj.calculate_d_range();

            % NUFFT for ISL and dISL/da
            a_range = [0, pi];

            obj.f = oNufft.Type3(obj.d_range, a_range, "epsilon",obj.params.tol/40);
            obj.df = oNufft.Type3(obj.d_range, a_range, "epsilon",obj.params.tol/40);
            obj.f.apply(obj.d, obj.g);
            obj.df.apply(obj.d, 1i*obj.d.*obj.g);
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

            if ~isnan(options.a_range)
                if obj.params.nr_candidates > 1
                    a_grid = linspace(options.a_range(1),options.a_range(2), obj.params.nr_candidates);
                else
                    a_grid = a(3);
                end
                
                a_range = options.a_range;
            else
                a_grid = 0:dA:pi;
                a_range = [a_grid(1), a_grid(end)];
            end
            
            % Select multiple MM candidates
            [~, idx_a] = min(real(obj.f.interpolate(a_grid)));
            a_mm = a_grid(idx_a);

            logger.log( ...
                "isl",  obj.f, ...
                "a",    a, ...
                "a_MM", a_mm, ...
                "i_MM", 0, ...
                "p",    3, ...
                "i_cd", 1 ...
                );

            % Majorize / Minimize
            for i = 1:obj.params.max_mm_iters
                tangent = real(obj.df.interpolate(a_mm)).';
                a_step = - tangent/obj.smoothness;
                
                if ~allbetween(a_mm + a_step, a_range(1), a_range(2))
                    break; % Out of reach for save NUFFT interpolation
                end

                if all([abs(a_step) < obj.params.tol, i > 1])
                    % fprintf("%d\n", i);
                    break; % No more progress
                end

                a_mm = a_mm + a_step;
                logger.log( ...
                    "a_MM", a_mm, ...
                    "i_MM", i, ...
                    "p",    3, ...
                    "i_cd", 1 ...
                    );
            end

            % Select final minima
            isl     = real(obj.f.interpolate(a_mm));
		    [~, k]  = min(isl);
            a(3)    = a_mm(k);

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

            % if obj.params.nr_candidates > 1
            %     [~, idx_a] = findpeaks(- real(f.interpolate(candidates)), "NPeaks", obj.params.nr_candidates, "SortStr","descend");
            %     if numel(idx_a) < obj.params.nr_candidates
            %         warning("CDMM: Could not find the required nr of MM candidates");
            %     end
            %     a_mm        = [ candidates(idx_a) ];
            % else
                [~, idx_a] = min(real(f.interpolate(candidates)));
                a_mm = candidates(idx_a);
            % end
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