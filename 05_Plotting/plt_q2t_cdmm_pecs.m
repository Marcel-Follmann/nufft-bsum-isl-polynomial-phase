close all; clearvars; clc; 

figure_handler = Figurehandler.Create("savePath","01_Figures\");
fig = {};

data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\q2t_cdmm_pecs_inf\");

% Omit Data, where time is larger than threshold and goal ISL was not
% reached
idx = data.isl_pecs <= data.isl_cdmm - 3; % Account for 3 dB Tolerance;
data.time_pecs(idx) = NaN;

mu = groupsummary(data, ["nr_samples", "poly_order"], @(x) mean(x, "all", "omitmissing"), ["time_cdmm", "time_pecs"]);
mu.Properties.VariableNames(end-1:end) = ["mean_time_cdmm", "mean_time_pecs"];

poly_order_vec = sort(unique(mu.poly_order), "descend") - 1;
nr_samples_vec = sort(unique(mu.nr_samples), "ascend");

groups = findgroups(data.poly_order);
n = sum(groups == 1);
success_rate = splitapply(@(x) {sum(isnan(x)) / n}, data.time_pecs, groups);
success_rate = [success_rate{:}].';

line_styles = [":", "--", "-"];
marker_styles = ["s", "^"];

fig{end+1} = figure("Name","Quality to Time");
for i = 1:numel(poly_order_vec)
    poly_order = poly_order_vec(i);
    h(1,i) = loglog( ...
        nr_samples_vec, ...
        mu.mean_time_cdmm(mu.poly_order == poly_order+1), ...
        line_styles(i) + marker_styles(1), ...
        "DisplayName", "BSUM $Q = "+num2str(poly_order) + "$"); hold on; 
    h(2,i) = loglog( ...
        nr_samples_vec, ...
        mu.mean_time_pecs(mu.poly_order == poly_order+1), ...
        line_styles(i) + marker_styles(2), ...
        "DisplayName", "PECS $Q = "+num2str(poly_order) + "$ $( " + num2str(round(success_rate(i),2)*100) + "\%) $"); hold on;
end

ylim([1e-4, 1e3]);
box on; %grid on; 

box on; % grid on
xticks(2.^(log2(nr_samples_vec(1)):log2(nr_samples_vec(end))))
lgd = legend(reshape(h.', [], 1), "Location","northwest", "NumColumns", 2, "EdgeColor","white");
lgd.set("IconColumnWidth", 15);
xlabel("Sequence Length (-)");
ylabel("Runtime (s)")
ax = gca;
ax.YGrid = 'on';
ax.YMinorGrid = 'off';

%% Formating
figure_handler.PrepareSvg(fig{1}, "colorList",["#D37611", "#005287"]);

figure_handler.GenerateOutputs()
clear ans;

