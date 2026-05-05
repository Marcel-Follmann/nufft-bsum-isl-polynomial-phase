close all; clearvars; clc; 

fh = Figurehandler.Create("savePath","01_Figures\");
f = {};
colorlist = flip(["#00934B", "#D37611", "#005287"]);

data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\rt_cdmm_fisl\");
idx = data.nr_samples == 2.^round(log2(data.nr_samples));
data = data(idx,:);
T = stack(data, ...
    ["time_cdmm_ood", "time_cdmm", "time_fisl"], ...
    "NewDataVariableName", "runtime", ...
    "IndexVariableName", "type");
nr_samples_vec = sort(unique(data.nr_samples), "ascend");

%% Visualization
f{end+1} = figure("Name","Runtime FISL vs CDMM A");
x = log2(T.nr_samples);
b = boxchart(x, T.runtime, 'GroupByColor',T.type, 'Notch','on');
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

lgd = legend(["BSUM-OOD", "BSUM", "FISL"], "Location","northwest", "Box","off");
lgd.IconColumnWidth = 10;

for i = 1:numel(b)
    b(i).JitterOutliers = 'on';
    b(i).MarkerStyle = '.';
    b(i).BoxFaceColor = hex2rgb(colorlist(i));
    b(i).MarkerColor = hex2rgb(colorlist(i));
end

fh.PreparePdf(f{1});
ax.YScale ="log";
fh.GenerateOutputs();