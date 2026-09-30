classdef FFT < handle

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
		function obj = FFT(params)
			arguments
				params Algorithm.Params
            end
			obj.params = params;

            [obj.d, obj.g] = calculate_gd(obj.params.nr_samples, obj.params.nr_params);


            obj.smoothness  = obj.calculate_smoothness();
            obj.d_range = obj.calculate_d_range();
            obj.periods = 2*pi./obj.calculate_gcd();
            obj.tol = obj.params.tol * sqrt(min(obj.smoothness)) / sqrt(obj.params.nr_params - 2);
        end

        function a = run(obj, a)
			arguments
				obj 
				a (:,1) double 
            end

            phase = a(3:end).'*obj.d;
            g_p = obj.g.*exp(1i*phase);
            gradient = Inf * ones(1, obj.params.nr_params);
            gradient(1:2) = 0;
            for i_cd = 1:obj.params.nr_cd_iterations
                for p = 3:obj.params.nr_params
                    % Bring into DFT format
                    d_p = obj.d(p-2,:);
                    phase = a(p).'*d_p;
                    g_p = g_p.*exp(-1i*phase);
                    dmax  = max(d_p);
                    ind   = d_p + dmax + 1;
                    n     = 2*dmax + 1;
                    d_fft = -dmax:dmax;
                    g_lin = accumarray(ind.', g_p.', [n 1]).';          % natürliche Reihenfolge -dmax..dmax
                    g_fft = ifftshift(g_lin .* exp(-1j*a(p)*d_fft));    % erst modulieren, dann shiften
                    f     = (-floor(n/2) : ceil(n/2)-1) / n;
                    
                    ISL = real(fftshift(fft(g_fft)));
                    [~, idx_a] = min(ISL);
                    a_mm = a(p) + 2*pi*f(idx_a);  
                     
                    % Majorize / Minimize
                    for i = 1:obj.params.max_mm_iters
                        
                        tangent = imag(sum(d_p .* g_p .* exp(1j*a_mm.*d_p)));
                        a_step = - tangent/obj.smoothness(p-2);
                        
                        gradient(p) = abs(tangent)  / sqrt(obj.smoothness(p-2)) ;
                        if gradient(p) < obj.tol
                            break; % No more progress
                        end

                        a_mm = a_mm + a_step;
                    end

                    % Select final minima
                    a(p)    = a_mm;

                    phase = a(p).'*d_p;
                    g_p = g_p.*exp(1i*phase);
                end

                if all(gradient) < obj.tol, break; end
            end
            %delete(logger);
        end

    end

    methods (Access = private)
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