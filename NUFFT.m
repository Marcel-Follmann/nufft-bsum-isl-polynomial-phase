classdef NUFFT < handle

    properties (Access = public)
        K (1,1)
        N (1,1)
        X (1,:)
        gamma (1,1)
        J (1,1)
        eta0 (1,1)
        Tinv
        C
        S
        f
        E
    end

    methods (Access = public)
        function obj = NUFFT(options)
            arguments
                options.x (:,1)
                options.K (1,1)
            end
            obj.K = options.K;
            obj.N = numel(options.x);
            
            obj.eta0 = (obj.N-1)/2;
            
            obj.E = sum(options.x.^2);

            pad = round((obj.K-obj.N)/2);
            x_padded = [zeros(pad, 1); options.x; zeros(pad, 1)];
            obj.K = numel(x_padded);
            obj.gamma = 2*pi/obj.K;

            obj.X = real(fft(ifftshift(x_padded)))/sqrt(obj.N);
            obj.f(:,1) = 0:1/obj.K:1-1/obj.K;
        end

        function e = error(obj, omega)
            arguments
                obj 
                omega (1,1,:) 
            end
            % Just compute the empirical worst case given in the paper
            e = zeros(numel(omega), 1);
            e(:,1) = sqrt(obj.N) * 0.75 * exp(-obj.J*( 0.29 + 1.03*log(obj.K/obj.N) ));
        end

        function prepareInterpolator(obj, options)
            arguments
                obj
                options.J
            end
            obj.J = options.J;
            obj.S = diag(ones(obj.N,1));
            n(:,1) = 0:(obj.N-1); 
            j(1,:) = 1:obj.J;
            obj.C = exp(1i*obj.gamma*(n-obj.eta0)*j)/sqrt(obj.N);

            kappa = j-j.';
            obj.Tinv = obj.Dirichlet(kappa);
        end

        function X_hat = interpolate(obj, omega)
            arguments
                obj 
                omega (1,1,:) 
            end
            indices = mod(obj.k0(omega) + (1:obj.J),obj.K) + 1;
            X_hat = pagemtimes(obj.X(indices),conj(obj.calculateInterpolators(omega)));
            X_hat = reshape(X_hat, numel(omega),1);
        end
    end

    properties (Access = private)
    end

    methods (Access = private)
        function km = k0(obj, omega)
            if mod(obj.J,2) == 0
                km = floor(omega/obj.gamma) - (obj.J)/2;
            else
                km = round(omega/obj.gamma) - (obj.J+1)/2;
            end
        end

        function u = calculateInterpolators(obj,omega)
            arguments
                obj 
                omega (1,1,:) 
            end
            j(:,1) = 1:obj.J;

            tmp = exp(-1i* ( omega - obj.gamma*(obj.k0(omega) + j) )*obj.eta0 );
            Lambda = zeros(obj.J, obj.J, numel(omega));
            idx = 1:obj.J+1:obj.J^2;
            ind = idx(:) + obj.J^2*(0:numel(omega)-1);
            Lambda(ind) = reshape(tmp,obj.J,numel(omega)); % Diagonalize



            kappa = omega/obj.gamma - obj.k0(omega) - j;
            r = obj.Dirichlet(kappa);

            u = pagemtimes( pagectranspose(Lambda) , pagemldivide( obj.Tinv , r) );
        end

        function y = Dirichlet(obj, kappa)
            mu = obj.K/obj.N;
            if obj.J >= 10-mu
                % fprintf("Use Large N approximation\n");
                y = sinc(kappa/mu);
            else
                y = 1/obj.N * sin(pi*kappa*obj.N/obj.K)./sin(pi*kappa/obj.K);
                t = kappa/obj.K;
                mask = abs(t - round(t)) < 10*eps(class(t));
                y(mask) = 1;
            end
        end
    end

end