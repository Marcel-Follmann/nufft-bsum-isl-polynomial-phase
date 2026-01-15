clc; close all; clear;

%% Initialization
N = 50;
n = 0:N-1;


%% Type of signal
% 'sinusoid' / 'quadratic'
sig_type = 'quadratic';
% Unwrapped Phase: psi & Wrapped phase: phi
switch sig_type
    case 'sinusoid'
        f0 = 1/N;
%         sf = 1e2;
        sf = 50;
        psi_u = sf*sin(2*pi*f0*n);
        x = exp(1j*psi_u);
    case 'quadratic'
%         a0 = -3; a1 = -0.5; a2 = -0.02;
%         a0 = 3; a1 = 0.5; a2 = 0.02;
%         a0 = 3; a1 = 5; a2 = 0.2;
%         a0 = 1; a1 = 1; a2 = 1;
%         a0 = 0; a1 = 0; a2 = 1;
%         a0 = 0; a1 = 1; a2 = 1;
%         a0 = 1; a1 = 0; a2 = 1;
        a0 = 1; a1 = 1; a2 = 1;

        psi_u = a0 + a1.*n + a2.*n.^2;
        x = exp(1j*psi_u);
end

phi_w = angle(x);
phi_c = wrap_phase(psi_u);
%% Visualization
figure(1); hold all; box on;
plot(psi_u);
plot(phi_w);
plot(phi_c, '--d');
legend('unwrapped phase','wrapped phase', 'wrapped cust')
xlabel('Sample Index');
ylabel('Phase [rad]');

%% Check the custom unwrap phase method

psi_c = unwrap_1d(phi_w);

psi_m = unwrap(phi_w);


figure(2); hold all; box on;
plot(psi_u, '-.o');
plot(psi_c, '-x');
plot(psi_m, '-.^');


xlabel('Sample Index');
ylabel('Phase [rad]');
legend('actual unwrapped phase','unwrapped phase cust', 'matlab unwrapped phase')


function unwrapped_phase = unwrap_1d(wrapped_phase)
% UNWRAP_1D Unwraps a 1D phase signal by removing 2*pi jumps
%   Input:  wrapped_phase (1D array in radians)
%   Output: unwrapped_phase (continuous phase)

    % Initialize output
    unwrapped_phase = wrapped_phase(:);  % Force column vector
    
    % Compute phase differences between adjacent samples
    dp = diff(unwrapped_phase);
    
    % Adjust jumps > pi by adding multiples of 2*pi
    dp_corrected = dp - 2*pi * round(dp / (2*pi));
    
    % Integrate corrections cumulatively
    cumulative_correction = cumsum(dp_corrected);
    
    % Apply correction to the original signal (skip first sample)
    unwrapped_phase(2:end) = unwrapped_phase(1) + cumulative_correction;
end

function wrapped_phase = wrap_phase(unwrapped_phase)
% WRAP_PHASE Converts an unwrapped phase signal back to the wrapped range [-π, π]
%   Input:  unwrapped_phase (1D array in radians)
%   Output: wrapped_phase (phase wrapped to [-π, π])

    % Force input to column vector (consistent handling)
    unwrapped_phase = unwrapped_phase(:);
    
    % Wrap phase using modulo-2π arithmetic
    wrapped_phase = mod(unwrapped_phase + pi, 2*pi) - pi;
    
    % Ensure exact π maps to π (not -π) for consistency
    wrapped_phase(abs(wrapped_phase - (-pi)) < eps(pi)) = pi;
end

function unwrapped_phase = my_unwrap(wrapped_phase)
% MY_UNWRAP Unwrap phase angle by detecting jumps greater than pi
%   unwrapped_phase = my_unwrap(wrapped_phase) unwraps the input phase vector/matrix
%   by detecting jumps greater than pi and correcting them with 2*pi offsets

% Initialize output with input phase
unwrapped_phase = wrapped_phase;
sz = size(wrapped_phase);

% Find the first non-singleton dimension for operation
dim = find(sz > 1, 1);
if isempty(dim)
    dim = 1;
end

% Calculate phase differences along the chosen dimension
dp = diff(unwrapped_phase, 1, dim);

% Calculate the phase difference modulo 2*pi
dp_mod = mod(dp + pi, 2*pi) - pi;

% Find where corrections are needed (jumps > pi)
dp_correct = dp_mod - dp;
dp_correct(abs(dp) < pi) = 0;

% Calculate cumulative corrections
cum_correction = cumsum(dp_correct, dim);

% Pad with zeros to match input size
correction = zeros(sz);
if dim == 1
    correction(2:end,:) = cum_correction;
elseif dim == 2
    correction(:,2:end) = cum_correction;
else
    % Handle higher dimensions if needed
    subs = repmat({':'}, 1, ndims(wrapped_phase));
    subs{dim} = 2:sz(dim);
    correction(subs{:}) = cum_correction;
end

% Apply correction
unwrapped_phase = unwrapped_phase + correction;
end