close all; clearvars; clc; 

rng("shuffle");

fh1 = Figurehandler.Create("savePath","01_Figures\");
fh2 = Figurehandler.Create("savePath","01_Figures\");
f = {};

params_2 = Algorithm.Params(...
    "nr_samples", 64, ...
    "nr_params",3 ...
    );

cdmm_2 = Algorithm.CDMM.OOD(params_2);
fisl = Algorithm.FISL(params_2);

[x0, a2] = params_2.get_random_sequence();

x_fisl = fisl.run(x0);
isl_fisl = params_2.compute_isl_dB(x_fisl);

a = cdmm_2.run(a2, "a_range",a2(3)+[-.5, .5]);
x_cdmm = params_2.generate_signal(a);
isl_cdmm = params_2.compute_isl_dB(x_cdmm);

fprintf("ISL FISL: %.2f dB\nISL CDMM: %.2f dB\n", isl_fisl, isl_cdmm);

r_fisl = params_2.compute_autocorr(x_fisl);
r_cdmm = params_2.compute_autocorr(x_cdmm);

N = params_2.nr_samples;
lags = (-(N-1):(N-1));
dNu = 1/N;
nuVec = -.5:dNu:.5-dNu;

AF_cdmm = pow2db(params_2.compute_ambiguity_function(x_cdmm, lags,nuVec));
AF_fisl = pow2db(params_2.compute_ambiguity_function(x_fisl, lags,nuVec));

params_3 = Algorithm.Params(...
    "nr_samples", 64, ...
    "nr_params",4 ...
    );
params_4 = Algorithm.Params(...
    "nr_samples", 64, ...
    "nr_params",5 ...
    );

cdmm_3 = Algorithm.CDMM.General(params_3);
cdmm_4 = Algorithm.CDMM.General(params_4);

MC = 10;
for mc = 1:MC
    [x0, a2] = params_2.get_random_sequence();
    [~, a3] = params_3.get_random_sequence();
    [~, a4] = params_4.get_random_sequence();
    x_fisl = fisl.run(x0);
    a_2 = cdmm_2.run(a2);
    a_3 = cdmm_3.run(a3);
    a_4 = cdmm_4.run(a4);
    x_cdmm_2 = params_2.generate_signal(a_2);
    x_cdmm_3 = params_3.generate_signal(a_3);
    x_cdmm_4 = params_4.generate_signal(a_4);

    Q_cdmm_2(mc) = pow2db(params_2.compute_doppler_tolerance(x_cdmm_2));
    Q_cdmm_3(mc) = pow2db(params_3.compute_doppler_tolerance(x_cdmm_3));
    Q_cdmm_4(mc) = pow2db(params_4.compute_doppler_tolerance(x_cdmm_4));
    Q_fisl(mc) = pow2db(params_4.compute_doppler_tolerance(x_fisl));
    fprintf(".");
end
fprintf("\nMean Q=2: %.2f dB\nMean Q=3: %.2f dB\nMean Q=4: %.2f dB\nMean FISL: %.2f dB\n", ...
    mean(Q_cdmm_2), mean(Q_cdmm_3), mean(Q_cdmm_4), mean(Q_fisl));
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
% hold on; plot(1/params.nr_samples * [1, 1], [-64, 64]);

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
