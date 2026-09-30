close all; clearvars; clc; 

fh = Figurehandler.Create("savePath","01_Figures\");
f = {};
colorlist = flip(["#00934B", "#D37611", "#005287"]);

data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\dt\");
% group = discretize(data.type, 1:5, "categorical", ["Q=2", "Q=3", "Q=4", "FISL"]);
% group = reordercats(group, {'Q = 2','Q = 3','Q = 4','FISL'});
f{end+1} = figure("Name","Doppler Tolerance");
b = boxchart(data.type, data.Q, "Notch","off");
ylabel("$\nu$ ( $dB$ )");
xticks(unique(data.type));
xticklabels(["Q=2", "Q=3", "Q=4", "FISL"])

box on;
for i = 1:numel(b)
    b(i).JitterOutliers = 'on';
    b(i).MarkerStyle = '.';
    b(i).BoxFaceColor = hex2rgb(colorlist(i));
    b(i).MarkerColor = hex2rgb(colorlist(i));
end

fh.PreparePdf(f{1});
fh.GenerateOutputs();

