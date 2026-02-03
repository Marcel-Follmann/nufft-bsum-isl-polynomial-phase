classdef PECS < handle

    properties
        nSamples
        nPolyOrder
        pNorm
        tol
        nMaxIters
        
        data_logger
    end

    properties (SetAccess = protected)
        
    end

    methods
        function obj = PECS(params)
            arguments
                params.nr_samples
                params.nr_params
                params.pNorm = 2
                params.tol = 1e-4
                params.max_mm_iters = 1e4
                
                params.save_data        (1,1) logical                   = false
                params.file_location    (1,1) string {mustBeFolder}     = "02_Data"
                params.file_name        (1,1) string                    = ""
            end

            obj.nSamples = params.nr_samples;
            obj.nPolyOrder = params.nr_params;
            obj.pNorm = params.pNorm;
            obj.tol = params.tol; 
            obj.nMaxIters = params.max_mm_iters;
            
            % Prepare data logger
            obj.data_logger = DataLogger("log_enabled",params.save_data, "file_location",params.file_location, "file_name",params.file_name);
            obj.data_logger.log("params", params);
        end

        function [xUpdate, aUpdate] = perfectSquare(obj, x, a)
            phaseReference = angle(x);
            magnitudeReference = abs(x);

            vSamples(:,1) = 1:obj.nSamples;
            V = vSamples .^(0:obj.nPolyOrder-1); % Vandermonde

            phaseModel = V*a;
            phaseError = phaseModel - phaseReference;

            % Taylor Subroutine
            wCos = (magnitudeReference/2)    .* cos(phaseError);
            wSin = magnitudeReference        .* sin(phaseError);
            b = -2*wCos.*phaseReference - 2*wCos.*phaseError + wSin;
            A = 2*wCos .* V;
            
            % Match Model
            [aUpdate, ~] = lsqr(A,-b, 1e-3);
            xUpdate = obj.generateSignal(aUpdate);
        end

        function [x0, a0] = getRandomSequence(obj)
            % Generate a sequence with a random code reference
            mCode = PhaseCodeType.Random;
            signal = mCode.generate(obj.nSamples);
            phaseReference = angle(signal);
            vSamples(:,1) = 1:obj.nSamples;

            polynomialModel = zeros(obj.nSamples, obj.nPolyOrder);
            for q = 1:obj.nPolyOrder
                polynomialModel(:,q) = vSamples.^(q-1);
            end
            [a0,~] = lsqr(polynomialModel, phaseReference);
            x0 = obj.generateSignal(a0);
        end

        function mpsl(obj, x0, a0)
            x = x0;
            a = a0;
            obj.data_logger.log("a", a);
            for i = 1:obj.nMaxIters
                xTmp = obj.mpslStep(x);
                [x, a] = obj.perfectSquare(xTmp, a);
                obj.data_logger.log("a", a);
            end
        end

        function delete(obj)
            delete(obj.data_logger);
        end

        function xUpdate = mpslStep(obj, x)
            arguments 
                obj,
                x (:,1)
            end
            F = fft(x, 2*obj.nSamples); % 2N for Aperiodic ACF
            rCirc = ifft(abs(F).^2);
            rPos = rCirc(2:obj.nSamples); % lags 1, .., nSamples-1
            rAbs = abs(rPos);

            % Normalize Sidelobes
            t = norm(rAbs, obj.pNorm) + eps;
            rho = rAbs ./ t;

            % Weights / Curvature terms
            a = (1 + (obj.nPolyOrder-1).*rho.^obj.nPolyOrder - obj.nPolyOrder.*rho.^(obj.nPolyOrder-1)) ./ (t - rAbs).^2;
            w = 0.5*obj.nPolyOrder/t^2 .* rho.^(obj.nPolyOrder-2);

            % Bound maximum Eigenvalue
            lagScale = (obj.nSamples-1:-1:1).';
            lambdaMaxL = max(lagScale .* a);

            wFlip = flipud(w);
            c = [0; w; 0; wFlip] .* rCirc; 
            mu = real(fft(c));

            % Bound Lambda u (odd/even split)
            lambdaU = 0.5*( max(mu(1:2:end)) + max(mu(2:2:end)) );

            % Apply operator via FFT
            y = ifft(mu .* F);

            % Update 
            alpha = obj.nSamples*lambdaMaxL + lambdaU;
            xUpdate = alpha*x - y(1:obj.nSamples);

        end

        function x = generateSignal(obj,a)
            arguments
                obj
                a (:, 1)
            end
            vSamples(:,1) = 1:obj.nSamples;
            V = vSamples .^(0:obj.nPolyOrder-1); % Vandermonde
            phase = V*a;
            x = exp(1i*phase);
        end

        function [isl_lin, isl_dB] = getIsl(obj)
            vSamples(:,1) = 1:obj.nSamples;
            V = vSamples .^(0:obj.nPolyOrder-1); % Vandermonde
            phase = V*obj.matObj.a;
        
	        x = exp(1i*phase);
	        
            F = fft(x, 2*obj.nSamples, 1); % 2N for Aperiodic ACF
            rCirc = ifft(abs(F).^2, [],1);
            rPos = rCirc(2:obj.nSamples,:); % lags 1, .., nSamples-1
            rAbs = abs(rPos)/obj.nSamples;

	        isl_lin = 2*sum(rAbs.^2,1)/obj.nSamples;
	        isl_dB = 10*log10(isl_lin);
        end
    end
end