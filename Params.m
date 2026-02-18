classdef Params
    properties (GetAccess = public, SetAccess = immutable)
        % Shared
        nr_samples
        nr_params
        max_mm_iters
        tol

        % CDMM specific
        nr_neighbors
        nr_candidates
        nr_cd_iterations
        oversampling

        % PECS specific
        p_norm
    end

    methods
        function obj = Params(params)
            arguments
				params.nr_samples       (1,1) double {mustBePositive}               % Sequence Length
				params.nr_params        (1,1) double {mustBePositive}               % Number of Parameters
                params.nr_neighbors     (1,1) double {mustBePositive}       = 4     % Number of neighbors to include for interpolation
                params.nr_candidates    (1,1) double {mustBePositive}       = 30    % Number of minima to pursue
                params.max_mm_iters     (1,1) double {mustBePositive}       = 1e3   % Maximum number of MM steps
                params.nr_cd_iterations (1,1) double {mustBePositive}       = 1     % Number of complete CD iterations
                params.oversampling     (1,1) double {mustBePositive}       = 2     % Oversampling for NUFFT
                params.tol              (1,1) double {mustBeNonnegative}    = 0     % Tolerance for local minima
                params.p_norm           (1,1) double {mustBeNonnegative}    = 2     % p-th Norm calculation
            end
            obj.nr_samples          = params.nr_samples;
			obj.nr_params           = params.nr_params;
            obj.nr_neighbors        = params.nr_neighbors;
            obj.max_mm_iters        = params.max_mm_iters;
            obj.nr_candidates       = params.nr_candidates;
            obj.nr_cd_iterations    = params.nr_cd_iterations;
            obj.oversampling        = params.oversampling;
            obj.tol                 = params.tol;
            obj.p_norm              = params.p_norm;
        end

        function [x0, a0] = get_random_sequence(obj)
            % Generate a sequence with a random code reference
            mCode = PhaseCodeType.Random;
            signal = mCode.generate(obj.nr_samples);
            phaseReference = unwrap(angle(signal));
            vSamples(:,1) = 1:obj.nr_samples;

            polynomialModel = zeros(obj.nr_samples, obj.nr_params);
            for q = 1:obj.nr_params
                polynomialModel(:,q) = vSamples.^(q-1);
            end
            [a0,~] = lsqr(polynomialModel, phaseReference);
            x0 = obj.generate_signal(a0);
        end

        function x = generate_signal(obj,a)
            arguments
                obj
                a (:, 1)
            end
            vSamples(:,1) = 1:obj.nr_samples;
            V = vSamples .^(0:obj.nr_params-1); % Vandermonde
            phase = V*a;
            x = exp(1i*phase);
        end

        function isl_dB = compute_isl_dB(obj,a)
			paramPhase = zeros(1, obj.nr_samples);
			for n=1:obj.nr_samples
				paramPhase(n) = 0;
				for q=1:obj.nr_params
					paramPhase(n) = paramPhase(n) + ...
						a(q)*n^(q-1);
				end
			end

			x = exp(1i*sum(paramPhase, 1));
			
			acfFull = xcorr(x)/obj.nr_samples;
			acf = acfFull(obj.nr_samples:end);
			isl_lin = 2*sum( abs(acf(2:end)).^2)/obj.nr_samples;
			isl_dB = 10*log10(isl_lin);
        end
    end
end