close all; clearvars; clc; 

addpath(".\");
addpath("01_Figures\");
addpath("02_Data\");
addpath("99_Utilities\utils\");

figure_handler = FigureHandler("savePath","01_Figures\");
fig = {};

cdmm = load("02_Data\cdmm_performance.mat");
pecs = load("02_Data\pecs_performance.mat");

% Reconstruct a over the Optimization
[a, a_mm] = reorder_cdmm(cdmm);
isl_dB_cdmm = get_isl(a, cdmm.params);
isl_dB_pecs = get_isl(pecs.a, pecs.params);

fig{end+1} = figure("Name", "ISL in Coordinate"); 
plot(cdmm.isl(1).f*2, 10*log10(cdmm.isl(1).X), '.', "DisplayName","ISL");hold on; 
plot(a_mm(1,:)/pi, isl_dB_cdmm(1,:), 'o', "DisplayName","Startposition");
plot(a_mm(:,1)/pi, isl_dB_cdmm(:,1), '.', "DisplayName","Path");
plot(a_mm(end,:)/pi, isl_dB_cdmm(end,:), 'x', "DisplayName","Endposition");
plot(a_mm(:,2:end)/pi, isl_dB_cdmm(:,2:end), '.', "HandleVisibility","off");
xlim([0,.5]);
xticks([0,.25, .5, .75,1, 1.25, 1.5, 1.75, 2]);
xticklabels({"$0$", "$\frac{\pi}{4}$", "$\frac{\pi}{2}$"})
xlabel("$a_p$ (-)");
ylabel("ISL (dB)");
legend("Location","northeast", "NumColumns", 2, "Box","off")
box on; %grid on; 

fig{end+1} = figure("Name", "Convergence Plot"); 
label = ["PECS: Initial", "CDMM: Initial", "CDMM: Candidate " + string(2:cdmm.params.nr_candidates)];
semilogx(isl_dB_pecs); hold on;
semilogx(isl_dB_cdmm);
box on; %grid on; 
ylabel("ISL (dB)");
xlabel("Number of MM steps (-)")
lgd = legend(label, "Location","best", "Box","off");
lgd.IconColumnWidth = 5;

fig{end+1} = figure("Name", "Convergence Plot Zoom"); 
label = "CDMM: Candidate " + string(2:cdmm.params.nr_candidates);
semilogx(isl_dB_cdmm(:,2:end));
box on; %grid on; 
ylabel("ISL (dB)");
xlabel("Number of MM steps (-)")
lgd = legend(label, "Location","best", "Box","off");
lgd.IconColumnWidth = 5;

figure_handler.PreparePdf(fig{1}, "colorList",["#005287", "#005287", "#005287", "#005287", "#D37611", "#005287", "#00934B", "#CCCCCC"])
for i = 2:numel(fig)
    figure_handler.PreparePdf(fig{i}, "colorList",flip(["#005287", "#D37611", "#00D16C", "#00934B", "#00753D", "#004725"]));
end
clear ans
figure_handler.GenerateOutputs();

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
    a = cdmm.a(:,1);
    a_mm = reshape(cdmm.a_MM, cdmm.params.nr_candidates, cdmm.params.max_mm_iters+1, cdmm.params.nr_params-2, cdmm.params.nr_cd_iterations);
    a_mm = permute(a_mm(:,:,1,1), [3, 2, 1]);
    a = a.*ones(cdmm.params.nr_params, cdmm.params.max_mm_iters+1, cdmm.params.nr_candidates);
    a(3,:,:) = a_mm;
    
    a_mm = permute(a_mm, [2, 3, 1]);
end