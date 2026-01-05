classdef NuFft < handle
    %NUFFT Summary of this class goes here
    %   Detailed explanation goes here

    properties
        K = NaN
        N = NaN
        J
        s 
        gamma
        T
        Lambda
        b
    end

    methods (Access = public)
        function obj = NuFft(options)
            arguments
                options.J
                options.s
            end
            obj.s = options.s;
            obj.J = options.J; 
        end

        function [Y, f] = kPointFft(obj, x, K)
            obj.N = numel(x);
            obj.K = K; 
            n(:,1) = -obj.N/2:obj.N/2-1; 
            
            Yk = zeros(obj.N, K/obj.N);
            for k = 0:K/obj.N-1
                Yk(:,k+1) = fft(fftshift(x.*exp(-1i*2*pi/K*n*k)));
            end
        
            Y = reshape(Yk.', [], 1);
            f(:,1) = -1/2:1/K:1/2-1/K; 
            
        end

        
        
    end

    methods (Access = public)
        function checkRequirements(obj)
            if isnan(obj.K)
                error("Run kPointFft(.) first");
            end
        end

        function prepareInterpolator(obj)
            obj.checkRequirements();
            obj.gamma = 2*pi/obj.K; 

            S = diag(obj.s);
            nu0 = (obj.N - 1) / 2;
            n(:,1) = 1:obj.N; 
            j(1,:) = 1:obj.J;
            C = exp(1i*obj.gamma*j.*(n-nu0))/sqrt(obj.N);

            obj.Lambda = @(omega) diag(exp(1i*(omega-obj.gamma*(obj.k0(omega) + j))*nu0));
            obj.b = @(omega) exp(1i*(omega - obj.gamma*obj.k0(omega)).*(n-nu0))/sqrt(obj.N);

            [Q, R] = qr(S'*C, "econ");
            obj.T = R\Q';
        end

        function u = interpolators(obj, omega)
            u = obj.Lambda(omega)' *obj.T * obj.b(omega);
        end

        function km = k0(obj,omega_m)
            obj.checkRequirements();

            
            
            if mod(obj.J,2) == 1 % odd
                km = round(omega_m/obj.gamma) - (obj.J+1)/2;
            else % even
                km = floor(omega_m/obj.gamma) - (obj.J)/2;
            end
        end
    end
end