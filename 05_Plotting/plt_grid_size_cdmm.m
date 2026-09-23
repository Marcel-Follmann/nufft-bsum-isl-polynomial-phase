close all; clearvars; clc;

fh = Figurehandler.Create("savePath","01_Figures\");
f = {};
colorlist = flip(["#00934B", "#D37611", "#005287"]);

folder_name = "grid_size";
data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\" + folder_name + "\");
load("02_Data\" + folder_name + "\settings.mat");

mu = groupsummary(data, ["nr_samples", "poly_order", "grid_size"], "mean", ["time", "isl"]);


[G, N, Q] = findgroups(data.nr_samples, data.poly_order); %categorical(data.poly_order);
valueset = unique(G, "sorted", "first");
catnames = "N = " + num2str(N) + ", Q = " + num2str(Q-1);
group = categorical(G, valueset, catnames);
% group = reordercats(group, {'Q = 2','Q = 3'});


f{end+1} = figure("Name", "ISL vs Grid Size");
x = log2(data.grid_size);
b = boxchart(x, data.isl, 'GroupByColor',group, 'Notch','off');
xt = log2(grid_size_vec);
xticks(xt);
xticklabels(string(grid_size_vec));
ds = (xt(2) - xt(1))/2;
xlim([xt(1) - ds, xt(end) + ds]);
ylim([-26, -12]);
lgd = legend("Location","northeast", "EdgeColor","white");
xlabel("Nr. of grid points ( $-$ )")
ylabel("ISL ( $dB$ )")
ax = gca;
ax.YGrid = 'on';
ax.YMinorGrid = 'off';
box on; 
lgd.IconColumnWidth = 10;
for i = 1:numel(b)
    b(i).JitterOutliers = 'on';
    b(i).MarkerStyle = '.';
    b(i).BoxFaceColor = hex2rgb(colorlist(i));
    b(i).MarkerColor = hex2rgb(colorlist(i));
end

f{end+1} = figure("Name", "Time vs Grid Size");
b = boxchart(x, data.time, 'GroupByColor',group, 'Notch','off');
xt = log2(grid_size_vec);
xticks(xt);
xticklabels(string(grid_size_vec));
ds = (xt(2) - xt(1))/2;
xlim([xt(1) - ds, xt(end) + ds]);
ylim([.5e-2, 1e1])
lgd = legend("Location","northeast", "EdgeColor","white", "NumColumns", 1);
xlabel("Nr. of grid points ( $-$ )")
ylabel("Runtime ( $s$ )")
ax = gca;
ax.YGrid = 'on';
ax.YMinorGrid = 'off';
box on; 
lgd.IconColumnWidth = 10;
for i = 1:numel(b)
    b(i).JitterOutliers = 'on';
    b(i).MarkerStyle = '.';
    b(i).BoxFaceColor = hex2rgb(colorlist(i));
    b(i).MarkerColor = hex2rgb(colorlist(i));
end

for i = 1:numel(f)
    fh.PreparePdf(f{i});
end
ax.YScale ="log";
fh.GenerateOutputs();