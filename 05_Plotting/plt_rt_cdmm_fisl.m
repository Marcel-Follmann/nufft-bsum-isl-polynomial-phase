close all; clearvars; clc; 

fh = Figurehandler.Create("savePath","01_Figures\");
f = {};
colorlist = flip(["#00934B", "#D37611", "#005287"]);
nr_files = 500;

data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\rt_cdmm_fisl\", nr_files);
group = findgroups(data.nr_samples); 

nr_monte_carlo = sum(group == 1);
mc = 1:nr_monte_carlo;

cummean_rt_cdmm = cummean(data.time_cdmm_ood, group);
cumstd_rt_cdmm = cumstd(data.time_cdmm_ood, group);

mean_convergence = abs(cummean_rt_cdmm - cummean_rt_cdmm(:, end));
convergence = max(mean_convergence, [], 1);
confidence = max(abs(1.96 * cumstd_rt_cdmm ./ sqrt(mc)), [], 1);

idx = data.nr_samples == 2.^round(log2(data.nr_samples));
data = data(idx,:);
T = stack(data, ...
    ["time_cdmm_ood", "time_cdmm", "time_fisl"], ...
    "NewDataVariableName", "runtime", ...
    "IndexVariableName", "type");
nr_samples_vec = sort(unique(data.nr_samples), "ascend");

%% Visualization
f{end+1} = figure("Name", "Convergence rt cdmm vs fisl");
semilogy(mc, mean_convergence.', "HandleVisibility","off"); hold on; 
semilogy(mc, confidence, "DisplayName", "$95\%$ Confidence"); hold on;
semilogy(mc, convergence, "DisplayName", "Max. Envelope"); hold on;
xlim([1, nr_monte_carlo])
xlabel("Monte Carlo Size (-)");
ylabel("Error to Final (s)");
lgd = legend("Location","southwest", "EdgeColor","white");
lgd.set("IconColumnWidth", 15);
% ylim([1e-1, 1e2])
ax = gca;
ax.YGrid = 'on';
ax.YMinorGrid = 'off';


f{end+1} = figure("Name","Runtime FISL vs CDMM A");
x = log2(T.nr_samples);
b = boxchart(x, T.runtime, 'GroupByColor',T.type, 'Notch','off');
xlabel("Sequence Length (-)");
ylabel("Runtime (s)");
xt = log2(nr_samples_vec);
xticks(xt);
xticklabels(string(nr_samples_vec));
ds = (xt(2) - xt(1))/2;
% xticks(nr_samples_vec)
xlim([xt(1) - ds, xt(end) + ds]);
ax = gca;
ax.YGrid = 'on';
ax.YMinorGrid = 'off';
box on; 

lgd = legend(["BSUM-OOD", "BSUM", "FISL"], "Location","northwest", "EdgeColor","white");
lgd.IconColumnWidth = 10;

for i = 1:numel(b)
    b(i).JitterOutliers = 'on';
    b(i).MarkerStyle = '.';
    b(i).BoxFaceColor = hex2rgb(colorlist(i));
    b(i).MarkerColor = hex2rgb(colorlist(i));
end

fh.PrepareSvg(f{1}, "colorList",["#005287", "#D37611", repelem("#CCCCCC", numel(mean_convergence))], "width", 31.56/2);
fh.PrepareSvg(f{2}, "width", 31.56/2);
ax.YScale ="log";
fh.GenerateOutputs();