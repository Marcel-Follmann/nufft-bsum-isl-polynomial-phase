close all; clearvars; clc; 

figure_handler = Figurehandler.Create("savePath","01_Figures\");
fig = {};

data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\t2q_cdmm_pecs\");
mu = groupsummary(data, ["nr_samples", "poly_order"], "mean", ["isl_cdmm", "isl_pecs"]);

poly_order_vec = sort(unique(mu.poly_order), "descend") - 1;
nr_samples_vec = sort(unique(mu.nr_samples), "ascend");

nr_monte_carlo = height(data) / height(mu);
mc = 1:nr_monte_carlo;

group = findgroups(data.nr_samples, data.poly_order);
cummean_isl_cdmm = cummean(data.isl_cdmm, group);
cumstd_isl_cdmm = cumstd(data.isl_cdmm, group);

mean_convergence = abs(cummean_isl_cdmm - cummean_isl_cdmm(:,end));
convergence = max(mean_convergence, [], 1);
confidence = max(abs(1.96*cumstd_isl_cdmm),[], 1);


fig{end+1} = figure("Name","Convergence Time 2 Quality");
semilogy(mc, mean_convergence, "HandleVisibility","off"); hold on; 
semilogy(mc, confidence, "DisplayName", "$95\%$ Confidence"); hold on;
semilogy(mc, convergence, "DisplayName", "Max. Envelope"); hold on;
xlim([1, nr_monte_carlo])
xlabel("Monte Carlo Size (-)");
ylabel("Error to Final (dB)");
lgd = legend("Location","northeast", "EdgeColor","white");
lgd.set("IconColumnWidth", 15);
ylim([1e-1, 1e2])
ax = gca;
ax.YGrid = 'on';
ax.YMinorGrid = 'off';
% ax.XGrid = 'on';
% ax.XMinorGrid = 'off';


line_styles = [":", "--", "-"];
marker_styles = ["s", "^"];

fig{end+1} = figure("Name","Time to Quality");
for i = 1:numel(poly_order_vec)
    poly_order = poly_order_vec(i);
    h(1,i) = semilogx( ...
        nr_samples_vec, ...
        mu.mean_isl_cdmm(mu.poly_order == poly_order+1), ...
        line_styles(i) + marker_styles(1), ...
        "DisplayName", "BSUM $Q = "+num2str(poly_order) + "$"); hold on; 
    h(2,i) = semilogx( ...
        nr_samples_vec, ...
        mu.mean_isl_pecs(mu.poly_order == poly_order+1), ...
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
figure_handler.PreparePdf(fig{1}, "colorList",["#005287", "#D37611", repelem("#CCCCCC", numel(mean_convergence))]);
figure_handler.PreparePdf(fig{2}, "colorList",["#D37611", "#005287"]);

figure_handler.GenerateOutputs()
clear ans;