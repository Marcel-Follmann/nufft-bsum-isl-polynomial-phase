close all; clearvars; clc; 

rng("shuffle");

fh1 = Figurehandler.Create("savePath","01_Figures\");
fh2 = Figurehandler.Create("savePath","01_Figures\");
f = {};

data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\rt_cdmm_fisl\");
mu = groupsummary(data, "nr_samples", "mean", ["time_cdmm_ood", "time_cdmm", "time_fisl"]);
nr_samples_vec = sort(unique(mu.nr_samples), "ascend");

params = Algorithm.Params(...
    "nr_samples", 64, ...
    "nr_params",3, ...
    "max_mm_iters",1e3, ...
    "tol", 1e-9, ...
    "nr_candidates",1000 ...
    );

cdmm = Algorithm.CDMM.OOD(params);
fisl = Algorithm.FISL(params);

[x0, a0] = params.get_random_sequence();
a0 = mod(a0, pi); % as this is periodic

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
nuVec = linspace(-.5, .5, 1e2+1);

AF_cdmm = 20*log10(params.compute_ambiguity_function(x_cdmm, lags,nuVec));
AF_fisl = 20*log10(params.compute_ambiguity_function(x_fisl, lags,nuVec));
% AF_cdmm(AF_cdmm <= -50) = NaN;
% AF_fisl(AF_fisl <= -50) = NaN;


% f{end+1} = figure("Name","AF");
% t = tiledlayout(1,2);
f{end+1} = figure("Name","AF BSUM");%nexttile;
imagesc(nuVec, lags, AF_cdmm);
clim([-20, 0]);
shading interp;
grid off;
map = CreateColormap((["#005287", "#FFFFFF", "#D37611"]));

ylabel("Lag (-)")
xlabel("Normalized Doppler (-)");
zlabel("Magnitude (dB)");
label_cdmm = sprintf("BSUM, ISL: %.2f dB", isl_cdmm);
% title(label_cdmm)
map = CreateColormap((["#005287", "#FFFFFF", "#D37611"]));
% colormap(map)
% c = colorbar;
clim([-20, 0]);
% ylabel(c, "Magnitude (dB)");
label_fisl = sprintf("FISL, ISL: %.2f dB", isl_fisl);
% title(label_fisl)
ax = gca;
ax.XMinorTick = 'off';
ax.YGrid ="on";

% f{end+1} = figure("Name","AF FISL");
f{end+1} = figure("Name","AF FISL");% nexttile;
imagesc(nuVec, lags, AF_fisl);
grid off;
% shading interp;
ylabel("Lag (-)")
xlabel("Normalized Doppler (-)");
zlabel("Magnitude (dB)");
map = CreateColormap((["#005287", "#FFFFFF", "#D37611"]));
% colormap(map)
% c = colorbar;
% c.Layout.Tile = "north";
clim([-20, 0]);
% ylabel(c, "Magnitude (dB)");
label_fisl = sprintf("FISL, ISL: %.2f dB", isl_fisl);
% title(label_fisl)
ax = gca;
ax.XMinorTick = 'off';
ax.YGrid ="on";



f{end+1} = figure("Name","Autocorr FISL vs CDMM");
plot(lags, 20*log10(abs(r_fisl)), '-', "DisplayName", label_fisl); hold on; 
plot(lags, 20*log10(abs(r_cdmm)), '-', "DisplayName", label_cdmm);
xlabel('Lag (-)');
ylabel('Magnitude (dB)');
grid off; 
lgd = legend("Location", "northeast", "EdgeColor","white");
clim([-20, 0])
lgd.IconColumnWidth = 10;
ax = gca;
ax.XMinorTick = 'off';
ax.YGrid ="on";

fh1.PrepareSvg(f{1}, "width", 31.56/2);
fh1.PrepareSvg(f{2}, "width", 31.56/2);
fh1.GenerateOutputs();
% fh2.GenerateOutputs();

function map = CreateColormap(colorlist)
    N = 256;

    for iColor = 1:numel(colorlist)
        rgbList(iColor, :) = hex2rgb(colorlist(iColor));
    end

    map = interp1(linspace(0,1,iColor), rgbList, linspace(0,1,N));
end
