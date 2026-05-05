close all; clearvars; clc; 

figure_handler = Figurehandler.Create("savePath","01_Figures\");
fig = {};

data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\rt_cdmm_pecs\");
mu = groupsummary(data, ["nr_samples", "poly_order"], "mean", ["timeCDMM", "timePECS"]);

poly_order_vec = sort(unique(mu.poly_order), "descend") - 1;
nr_samples_vec = sort(unique(mu.nr_samples), "ascend");

line_styles = [":", "--", "-"];
marker_styles = ["s", "^"];

fig{end+1} = figure("Name","Runtime Samples");
for i = 1:numel(poly_order_vec)
    poly_order = poly_order_vec(i);
    h(1,i) = semilogx( ...
        nr_samples_vec, ...
        mu.mean_timeCDMM(mu.poly_order == poly_order+1), ...
        line_styles(i) + marker_styles(1), ...
        "DisplayName", "BSUM $Q = "+num2str(poly_order) + "$"); hold on; 
    h(2,i) = semilogx( ...
        nr_samples_vec, ...
        mu.mean_timePECS(mu.poly_order == poly_order+1), ...
        line_styles(i) + marker_styles(2), ...
        "DisplayName", "PECS $Q = "+num2str(poly_order) + "$"); hold on;
end
% ylim([-35, 0]);
box on; % grid on
xticks(2.^(log2(nr_samples_vec(1)):log2(nr_samples_vec(end))))
lgd = legend(reshape(h.', [], 1), "Location","northeast", "NumColumns", 2, "EdgeColor","white");
lgd.set("IconColumnWidth", 15);
xlabel("Sequence Length (-)");
ylabel("ISL (dB)")
ax = gca;
ax.YGrid = 'on';
ax.YMinorGrid = 'off';
% ax.XGrid = 'on';
% ax.XMinorGrid = 'off';


%% Formating
figure_handler.PreparePdf(fig{1}, "colorList",["#D37611", "#005287"]);

figure_handler.GenerateOutputs()
clear ans;