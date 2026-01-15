classdef NUFFT < handle

    properties (Access = public)
        K (1,1)
        N (1,1)
        X (1,:)
        gamma (1,1)
        J (1,1)
        eta0 (1,1)
        T
        C
        S
    end

    methods (Access = public)
        function obj = NUFFT(options)
            arguments
                options.x (:,1)
                options.K (1,1)
            end
            obj.K = options.K;
            obj.N = numel(options.x);
            obj.gamma = 2*pi/obj.K;
            obj.eta0 = (obj.N-1)/2;
            obj.X = fftshift(ifft(options.x, obj.K))*obj.K;
        end

        function e = error(obj, omega)
            j(:,1) = 1:obj.J;
            r = sinc(omega/obj.gamma - obj.k0(omega) - j);
            e = sqrt(1 - r.'*obj.T*r)*sqrt(obj.N);
        end

        function prepareInterpolator(obj, options)
            arguments
                obj
                options.J
            end
            obj.J = options.J;
            obj.S = diag(ones(obj.N,1));
            n(:,1) = 1:obj.N; 
            j(1,:) = 1:obj.J;
            obj.C = exp(1i*obj.gamma*(n-obj.eta0)*j)/sqrt(obj.N);

            tmp = sinc(j-j.');
            obj.T = inv(tmp);
        end

        function X_hat = interpolate(obj, omega)
            indices(1,:) = mod(obj.k0(omega) + (1:obj.J),obj.K) + 1;
            X_hat = obj.X(indices)*conj(obj.calculateInterpolators(omega));
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
            j(:,1) = 1:obj.J;
            Lambda = diag([exp(-1i* ( omega - obj.gamma*(obj.k0(omega) + j) )*obj.eta0 )]);
            mu = obj.K/obj.N;

            r = sinc(omega/obj.gamma - obj.k0(omega) - j);

            u = Lambda*obj.T*r;
        end
    end

end