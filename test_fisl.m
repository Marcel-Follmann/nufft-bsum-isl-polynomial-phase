close all; clearvars; clc; 

addpath(".\");
addpath("02_Data\");
addpath("99_Utilities");
rng("shuffle");

fh = Figurehandler.Create("savePath","01_Figures\");
f = {};

params = Params(...
    "nr_samples", 128, ...
    "nr_params",3, ...
    "max_mm_iters",1e3, ...
    "tol", 1e-5, ...
    "nr_candidates",5 ...
    );

cdmm = CDMM.Quadratic(params);
fisl = FISL(params);

[x0, a0] = params.get_random_sequence();
a0 = mod(a0, 2*pi) + pi; % as this is periodic

tic; 
x_fisl = fisl.run(x0);
time_fisl = toc;
fprintf("Time FISL: %.2e s\n", time_fisl);
isl_fisl = params.compute_isl_dB(x_fisl);

tic; 
a = cdmm.run(a0, "a_range",a0(3)+[-.5, .5]);
time_cdmm = toc; 
fprintf("Time CDMM: %.2e s\n", time_cdmm);

x_cdmm = params.generate_signal(a);
isl_cdmm = params.compute_isl_dB(x_cdmm);

fprintf("ISL FISL: %.2f dB\nISL CDMM: %.2f dB\n", isl_fisl, isl_cdmm);

r_fisl = params.compute_autocorr(x_fisl);
r_cdmm = params.compute_autocorr(x_cdmm);


lags = (-(params.nr_samples-1):(params.nr_samples-1));
nuVec = linspace(-.5, .5, 1e3+1);

AF_cdmm = 20*log10(params.compute_ambiguity_function(x_cdmm, lags,nuVec));
AF_fisl = 20*log10(params.compute_ambiguity_function(x_fisl, lags,nuVec));
AF_cdmm(AF_cdmm <= -50) = NaN;
AF_fisl(AF_fisl <= -50) = NaN;


f{end+1} = figure("Name","AF CDMM");
surf(nuVec, lags, AF_cdmm);
zlim([-50, 0]);
clim([-30, 0]);
shading interp;
ylabel("Lag (-)")
xlabel("Normalized Doppler (-)");
zlabel("Magnitude (dB)");

f{end+1} = figure("Name","AF FISL");
surf(nuVec, lags, AF_fisl);
zlim([-50, 0]);
clim([-30, 0]);
shading interp;
ylabel("Lag (-)")
xlabel("Normalized Doppler (-)");
zlabel("Magnitude (dB)");



f{end+1} = figure("Name","Autocorr FISL vs CDMM");
plot(lags, 20*log10(abs(r_fisl)), '-', "DisplayName", "FISL"); hold on; 
plot(lags, 20*log10(abs(r_cdmm)), '-', "DisplayName", "CDMM");
xlabel('Lag');
ylabel('Aperiodic ACF sidelobes (dB)');
grid on;
legend("Location", "best");

for i = 1:numel(f)
    fh.PrepareSvg(f{i});
end
% fh.GenerateOutputs();

