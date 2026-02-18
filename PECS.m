classdef PECS < handle

    properties
        params Params = Params.empty
    end

    properties (SetAccess = protected)
        
    end

    methods
        function obj = PECS(params)
            arguments
                params Params = Params.empty
            end
            obj.params = params;
        end

        function [xUpdate, aUpdate] = perfectSquare(obj, x, a)
            phaseReference = angle(x);
            magnitudeReference = abs(x);

            vSamples(:,1) = 1:obj.params.nr_samples;
            V = vSamples .^(0:obj.params.nr_params-1); % Vandermonde

            phaseModel = V*a;
            phaseError = phaseModel - phaseReference;

            % Taylor Subroutine
            wCos = (magnitudeReference/2)    .* cos(phaseError);
            wSin = magnitudeReference        .* sin(phaseError);
            b = -2*wCos.*phaseReference - 2*wCos.*phaseError + wSin;
            A = 2*wCos .* V;
            
            % Match Model
            [aUpdate, flag] = lsqr(A,-b, 1e-3);
            if flag > 0
                msg = sprintf("Tolerance for lsqr calculation too tight. Error Code: %d", flag);
                warning("PECS:Perfect_Square_Calculation", msg);
            end
            xUpdate = obj.params.generate_signal(aUpdate);
        end

        function mpsl(obj, x0, a0, options)
            arguments
                obj, 
                x0, 
                a0, 
                options.data_logger Datalogger.Type = Datalogger.Empty()
            end
            options.data_logger.log("params", obj.params);

            x = x0;
            a = a0;
            options.data_logger.log("a", a);
            for i = 1:obj.params.max_mm_iters
                xTmp = obj.mpslStep(x);
                [x, a_] = obj.perfectSquare(xTmp, a);
                if norm(a_ - a) < obj.params.tol
                    break;
                end
                a = a_;
                options.data_logger.log("a", a);
            end
            delete(options.data_logger);
        end

        function xUpdate = mpslStep(obj, x)
            arguments 
                obj,
                x (:,1)
            end
            F = fft(x, 2*obj.params.nr_samples); % 2N for Aperiodic ACF
            rCirc = ifft(abs(F).^2);
            rPos = rCirc(2:obj.params.nr_samples); % lags 1, .., nr_samples-1
            rAbs = abs(rPos);

            % Normalize Sidelobes
            t = norm(rAbs, obj.params.p_norm) + eps;
            rho = rAbs ./ t;

            % Weights / Curvature terms
            a = (1 + (obj.params.nr_params-1).*rho.^obj.params.nr_params - obj.params.nr_params.*rho.^(obj.params.nr_params-1)) ./ (t - rAbs).^2;
            w = 0.5*obj.params.nr_params/t^2 .* rho.^(obj.params.nr_params-2);

            % Bound maximum Eigenvalue
            lagScale = (obj.params.nr_samples-1:-1:1).';
            lambdaMaxL = max(lagScale .* a);

            wFlip = flipud(w);
            c = [0; w; 0; wFlip] .* rCirc; 
            mu = real(fft(c));

            % Bound Lambda u (odd/even split)
            lambdaU = 0.5*( max(mu(1:2:end)) + max(mu(2:2:end)) );

            % Apply operator via FFT
            y = ifft(mu .* F);

            % Update 
            alpha = obj.params.nr_samples*lambdaMaxL + lambdaU;
            xUpdate = alpha*x - y(1:obj.params.nr_samples);
        end
    end
end