close all; clearvars; clc; 

addpath(".\");
addpath("02_Data\");
addpath("99_Utilities");
rng("shuffle");






P = 32;
opts.maxIter = 1e3;
opts.tol = 1e-5;
opts.verbose = false;

n_poly_order = 3;

params = Params(...
    "nr_samples", P, ...
    "nr_params",n_poly_order, ...
    "max_mm_iters",opts.maxIter, ...
    "tol", opts.tol, ...
    "nr_candidates",3 ...
    );

cdmm = CDMM.Quadratic(params);

[x0, a0] = params.get_random_sequence();
a0 = mod(a0, 2*pi) + pi; % as this is periodic

tic; 
[z, hist] = fisl(P, x0, opts);
time_fisl = toc;
fprintf("Time FISL: %.2e s\n", time_fisl);

tic; 
a = cdmm.run(a0, "a_range",a0(3)+[-.5, .5]);
time_cdmm = toc; 
fprintf("Time CDMM: %.2e s\n", time_cdmm);

isl_cdmm = params.compute_isl_dB(a);
figure;
plot(10*log10(hist.ISL), 'LineWidth', 1.5); hold on; 
plot([1,numel(hist.ISL)], [isl_cdmm,isl_cdmm])
xlabel('Iteration');
ylabel('ISL (dB)');
grid on;


x = params.generate_signal(a);
% 
% 
% X = x .* D;                           % jede Spalte: x mit anderem Doppler
% Z = z .* D;
% 
% 

r1 = aperiodic_autocorr_from_fft(z);
r2 = aperiodic_autocorr_from_fft(x);


tauVec = (-(P-1):(P-1));
nuVec = linspace(-.5, .5, 1e3+1);

AFx = ambiguity_function(x, tauVec,nuVec);
AFz = ambiguity_function(z, tauVec,nuVec);
figure();
surf(nuVec, tauVec, 20*log10(AFx));
% view(0,90);
zlim([-50, 0]);
clim([-50, 0]);
% clim(zlim);
shading interp;
ylabel("Lag (-)")
xlabel("Normalized Doppler (-)");
zlabel("Magnitude (dB)");

figure();
surf(nuVec, tauVec, 20*log10(AFz));
% view(0,90);
zlim([-50, 0]);
clim([-50, 0]);
% clim(zlim);
shading interp;
ylabel("Lag (-)")
xlabel("Normalized Doppler (-)");
zlabel("Magnitude (dB)");


% 
% RD1 = fftshift(fft(r1, 1000, 2), 2);
% RD2 = fftshift(fft(r2, 1000, 2), 2);
% 
lags = 0:P-1;
figure;
plot(lags, 20*log10(abs(r1)), '-'); hold on; 
plot(lags, 20*log10(abs(r2)), '-');
xlabel('Lag');
ylabel('Aperiodic ACF sidelobes (dB)');
grid on;

function r = aperiodic_autocorr_from_fft(z)
% Returns r = [r(0); r(1); ...; r(P-1)]
    P = length(z);
    zp = [z; zeros(P,1)];
    c = ifft(abs(fft(zp)).^2);
    r = c(1:P)/P;
end

function AFn = ambiguity_function(x, tauVec, nuVec)
    P = numel(x);
    AF = zeros(length(tauVec), length(nuVec));
    n = (0:P-1).';
    for itau = 1:length(tauVec)
        tau = tauVec(itau);
        
        if tau >= 0
            x1 = x(1+tau:P);
            x2 = x(1:P-tau);
            n_use = n(1:length(x1));
        else
            t = -tau;
            x1 = x(1:P-t);
            x2 = x(1+t:P);
            n_use = n(1:length(x1));
        end
        
        for inu = 1:length(nuVec)
            nu = nuVec(inu);
            AF(itau,inu) = sum(x1 .* conj(x2) .* exp(1j*2*pi*nu*n_use));
        end
    end
    AFn = abs(AF) / max(abs(AF(:)));
end