classdef FISL < handle
% FISL for unimodular sequence design by ISL minimization
% Reference:
% Sankuru, Babu, "A fast iterative algorithm to design phase-only sequences
% by minimizing the ISL metric" (FISL). Uses BEFFT bound and FFT/IFFT
% implementation. See Algorithm 4 and Eqs. (25)-(46).

    properties (Access = private)
        params Algorithm.Params
    end

    methods (Access=public)
        function obj = FISL(params)
            arguments
				params Algorithm.Params
            end
            obj.params = params;
        end

        function z = run(obj, z0, options)
            arguments
                obj 
                z0 
                options.data_logger Datalogger.Type = Datalogger.Empty
            end
            logger = options.data_logger;
            
            z = exp(1j * angle(z0)); % enforce unimodular
            ISL_new = obj.params.compute_isl_dB(z);

            logger.log("params", obj.params, "z", z, "ISL", ISL_new);

            for iter = 1:obj.params.max_mm_iters
                r = obj.params.compute_autocorr(z, "symmetric",false);
                M_scalar = obj.befft_bound_from(r);
                yz = obj.toeplitz_R_times_z(r, z);

                % MM Update
                a_tilde = -(yz - M_scalar*z/4);

                z = exp(1j * angle(a_tilde)); % Project to unit circle

                ISL_old = ISL_new;
                ISL_new = obj.params.compute_isl_dB(z);

                logger.log("z", z, "ISL", ISL_new);

                rel_change = abs(ISL_new - ISL_old) / max(1, ISL_old);
                if rel_change <= obj.params.tol
                    break;
                end
            end
        end
    end

    methods (Access=private)
        function Mscalar = befft_bound_from(obj, r)
        % Implements Eq. like M = 4*(max even bins + max odd bins)
        % where s is DFT of d = [r(0), r(1), ..., r(P-1), 0, r*(P-1), ..., r*(1)]^T
        
            d = [r; 0; conj(r(end:-1:2))];   % length 2P
            s = fft(d);                      % same scaling convention as paper up to FFT convention
        
            split_even_odd = reshape(s,2,[]);
            max_values = max(split_even_odd, [], 2);

            Mscalar = 4 * sum(max_values);
        end
        
        function yz = toeplitz_R_times_z(obj, r, z)
        % Computes yz = R(z) * z
        % R is Hermitian Toeplitz with first column [r(0); r(1); ...; r(P-1)]
        % and first row [r(0), conj(r(1)), ..., conj(r(P-1))]
        
            P = length(z);
            c = [r; 0; conj(r(end:-1:2))];   % circulant embedding, length 2P
            zp = [z; zeros(P,1)];
            yz_full = ifft(fft(c) .* fft(zp));
            yz = yz_full(1:P);
        end 
    end
end