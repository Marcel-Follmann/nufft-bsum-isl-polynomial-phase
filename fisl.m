function [z, hist] = fisl(P, z0, opts)
% FISL for unimodular sequence design by ISL minimization
%
% Inputs
%   P   : sequence length
%   z0  : initial unimodular sequence (Px1), optional
%   opts.maxIter : maximum iterations, default 2000
%   opts.tol     : relative ISL tolerance, default 1e-8
%   opts.verbose : true/false, default false
%
% Outputs
%   z    : designed unimodular sequence
%   hist : struct with ISL history and relative change
%
% Reference:
% Sankuru, Babu, "A fast iterative algorithm to design phase-only sequences
% by minimizing the ISL metric" (FISL). Uses BEFFT bound and FFT/IFFT
% implementation. See Algorithm 4 and Eqs. (25)-(46).

    if nargin < 2 || isempty(z0)
        z0 = exp(1j * 2*pi * rand(P,1));
    end
    if nargin < 3
        opts = struct();
    end
    if ~isfield(opts, 'maxIter'),  opts.maxIter = 2000;     end
    if ~isfield(opts, 'tol'),      opts.tol = 1e-8;         end
    if ~isfield(opts, 'verbose'),  opts.verbose = false;    end

    z = exp(1j * angle(z0(:)));   % enforce unimodular
    hist.ISL = zeros(opts.maxIter,1);
    hist.relchg = zeros(opts.maxIter,1);

    for k = 1:opts.maxIter
        % Step 1: aperiodic autocorrelation r(l), l=0,...,P-1
        r = aperiodic_autocorr_from_fft(z);   % length P, r(1)=r(0)

        % Step 2: BEFFT bound for M
        Mscalar = befft_bound_from_r(r, P);   % scalar M, i.e. M = Mscalar * I

        % Step 3: compute y = R(z) * z without explicitly forming R
        yz = toeplitz_R_times_z(r, z);

        % Step 4: majorization update
        a_tilde = -(yz - 0.25 * Mscalar * z);

        % Step 5: project to unit circle elementwise
        z_new = exp(1j * angle(a_tilde));

        % Step 6: evaluate ISL
        r_new = aperiodic_autocorr_from_fft(z_new);
        isl_new = sum(abs(r_new(2:end)).^2)/P;

        hist.ISL(k) = isl_new;
        if k > 1
            hist.relchg(k) = abs(hist.ISL(k) - hist.ISL(k-1)) / max(1, hist.ISL(k-1));
            if opts.verbose
                fprintf('Iter %4d: ISL = %.6e, relchg = %.3e\n', k, isl_new, hist.relchg(k));
            end
            if hist.relchg(k) <= opts.tol
                z = z_new;
                hist.ISL = hist.ISL(1:k);
                hist.relchg = hist.relchg(1:k);
                return;
            end
        else
            if opts.verbose
                fprintf('Iter %4d: ISL = %.6e\n', k, isl_new);
            end
        end

        z = z_new;
    end

    % trim if maxIter reached
    nz = find(hist.ISL ~= 0, 1, 'last');
    if ~isempty(nz)
        hist.ISL = hist.ISL(1:nz);
        hist.relchg = hist.relchg(1:nz);
    else
        hist.ISL = [];
        hist.relchg = [];
    end
end

function r = aperiodic_autocorr_from_fft(z)
% Returns r = [r(0); r(1); ...; r(P-1)]
    P = length(z);
    zp = [z; zeros(P,1)];
    c = ifft(abs(fft(zp)).^2);
    r = c(1:P)/P;
end

function Mscalar = befft_bound_from_r(r, P)
% Implements Eq. like M = 4*(max even bins + max odd bins)
% where s is DFT of d = [r(0), r(1), ..., r(P-1), 0, r*(P-1), ..., r*(1)]^T

    d = [r; 0; conj(r(end:-1:2))];   % length 2P
    s = fft(d);                      % same scaling convention as paper up to FFT convention

    odd_idx  = 1:2:(2*P);            % MATLAB 1-based
    even_idx = 2:2:(2*P);

    % s should be real up to numerical noise for Hermitian-symmetric d
    s_odd_max  = max(real(s(odd_idx)));
    s_even_max = max(real(s(even_idx)));

    Mscalar = 4 * (s_odd_max + s_even_max);
end

function yz = toeplitz_R_times_z(r, z)
% Computes yz = R(z) * z
% R is Hermitian Toeplitz with first column [r(0); r(1); ...; r(P-1)]
% and first row [r(0), conj(r(1)), ..., conj(r(P-1))]

    P = length(z);
    c = [r; 0; conj(r(end:-1:2))];   % circulant embedding, length 2P
    zp = [z; zeros(P,1)];
    yz_full = ifft(fft(c) .* fft(zp));
    yz = yz_full(1:P);
end