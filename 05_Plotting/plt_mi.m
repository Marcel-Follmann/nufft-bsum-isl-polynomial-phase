close all; clearvars; clc; 

fh = Figurehandler.Create("savePath","01_Figures\");
f = {};
colorlist = flip(["#00934B", "#D37611", "#005287"]);

loc = "mi_uo\";
data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\" + loc);
load("02_Data\" + loc + "settings.mat");

f{end+1} = figure("Name","Mutual Interference");
group = categorical(data.nr_samples);
group = reordercats(group, {'N = 4','N = 16','N = 64'});
b = boxchart(data.poly_order, 20*log10(data.mi), 'GroupByColor',group, 'Notch','on');
ylabel("Normed Max. MI (dB)");
xlabel("Polynomial Order (-)");
box on;
lgd = legend("Location","southeast", "EdgeColor","white", "NumColumns",3);
lgd.IconColumnWidth = 10;
for i = 1:numel(b)
    b(i).JitterOutliers = 'on';
    b(i).MarkerStyle = '.';
    b(i).BoxFaceColor = hex2rgb(colorlist(i));
    b(i).MarkerColor = hex2rgb(colorlist(i));
end
ax = gca;
ylim([-25, 0])
ds = (poly_order_vec(2) - poly_order_vec(1))/2;
xticks(poly_order_vec-1)
xlim([poly_order_vec(1) - ds, poly_order_vec(end) + ds]-1);

for i = 1:numel(f)
    fh.PreparePdf(f{i});
end
ax.XMinorTick = 'off';
ax.YGrid ="on";
fh.GenerateOutputs()