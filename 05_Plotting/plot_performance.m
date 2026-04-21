close all; clearvars; clc; 

figure_handler = Figurehandler.Create("savePath","01_Figures\");
fig = {};

cdmm = load("02_Data\cdmm_performance.mat");
pecs = load("02_Data\pecs_performance.mat");

% Reconstruct a over the Optimization
[a, a_mm] = reorder_cdmm(cdmm);
isl_dB_cdmm = get_isl(a, cdmm.params);
isl_dB_pecs = get_isl(pecs.a, pecs.params);

% a_mm = mod(a_mm, pi);

fig{end+1} = figure("Name", "ISL in Coordinate"); 
% theta = linspace(0, 2*pi, 1e4);
p = numel(cdmm.isl)+1;
theta = linspace(cdmm.isl(end).s_range(1), cdmm.isl(end).s_range(2), 1e6);
plot(a_mm(end,:)/pi, isl_dB_cdmm(end,:), 'x', "DisplayName","Minima"); hold on; 
plot(theta/(pi), 10*log10(real(cdmm.isl(end).interpolate(theta))), '-', "DisplayName","ISL");hold on; 
xlabel("$a_" + num2str(p) + " / \pi$ (-)");
ylabel("ISL (dB)");
lgd = legend("Location","northeast", "Box", "off", "NumColumns", 2);
lgd.IconColumnWidth = 10;
box on; %grid on; 

fig{end+1} = figure("Name", "Convergence Plot"); 
semilogx(isl_dB_pecs, '.', "DisplayName","PECS"); hold on;
semilogx(isl_dB_cdmm, '.', "DisplayName","CDMM");
box on; %grid on; 
ylabel("ISL (dB)");
xlabel("Number of MM steps (-)")
lgd = legend("Location","northeast", "Box","off", "NumColumns", 1);
lgd.IconColumnWidth = 10;

fig{end+1} = figure("Name","ACF");
N = cdmm.params.nr_samples;
f = -N:N-1;
r_start = get_signal(cdmm.a(:,1), cdmm.params);
r_cdmm = get_signal(cdmm.a(:,end), cdmm.params);
r_pecs = get_signal(pecs.a(:,end), pecs.params);
plot(f,20*log10(r_start), "DisplayName","Initial"); hold on; 
plot(f,20*log10(r_cdmm), "DisplayName","CDMM"); 
plot(f,20*log10(r_pecs), "DisplayName","PECS");
ylim([-65, 5]);
ylabel("Magnitude (dB)");
xlabel("Lags (-)");
lgd = legend("Location","northeast", "Box", "off", "NumColumns", 1);
lgd.IconColumnWidth = 10;

figure_handler.PreparePdf(fig{1});%, "colorList",["#005287", "#005287", "#005287", "#005287", "#D37611", "#005287", "#00934B", "#CCCCCC"])
for i = 2:numel(fig)
    figure_handler.PreparePdf(fig{i});%, "colorList",flip(["#005287", "#D37611", "#00D16C", "#00934B", "#00753D", "#004725"]));
end
clear ans
figure_handler.GenerateOutputs();

function rAbs = get_signal(a, params)
    samples(:,1) = 1:params.nr_samples;
    V = samples .^ (0:params.nr_params-1); % Vandermonde
    phase = pagemtimes(V,a);
    
    x = exp(1i*phase);
    F = fft(x, 2*params.nr_samples, 1); % 2N for Aperiodic ACF
    rCirc = ifftshift(ifft(abs(F).^2, [],1));
    rAbs = abs(rCirc)/params.nr_samples;
end

function isl_dB = get_isl(a, params)
    samples(:,1) = 1:params.nr_samples;
    V = samples .^ (0:params.nr_params-1); % Vandermonde
    phase = pagemtimes(V,a);
    
    x = exp(1i*phase);
    
    F = fft(x, 2*params.nr_samples, 1); % 2N for Aperiodic ACF
    rCirc = ifft(abs(F).^2, [],1);
    rPos = rCirc(2:params.nr_samples,:,:); % lags 1, .., nSamples-1
    rAbs = abs(rPos)/params.nr_samples;
    
    isl_lin = 2*sum(rAbs.^2,1)/params.nr_samples;
    isl_dB = 10*log10(isl_lin);
    
    isl_dB = permute(isl_dB, [2,3,1]);
end

function [a, a_mm] = reorder_cdmm(cdmm)
    a = cdmm.a(:,end-1);
    a_mm = cdmm.a_MM(:, find(cdmm.i_MM == 0,1,"last"):end);
    % a_mm = a_mm(:,:,end,end);
    a_mm = permute(a_mm, [3, 2, 1]);
    a = a.*ones(cdmm.params.nr_params, size(a_mm,2), 1, cdmm.params.nr_cd_iterations);
    a(end,:,:) = a_mm;
    
    a_mm = permute(a_mm, [2, 3, 1]);
end