close all; clearvars; clc; 

figure_handler = Figurehandler.Create("savePath","01_Figures\");
fig = {};

file = uigetfile('*.mat', "Select Data", "02_Data");
load("02_Data\" + string(file));

reorder = @(data) reorder_data(data, [3, 2, 1], ...
    nr_monte_carlos, poly_order, samples); % Resembles for loop order in experiments file

cdmm = reorder(cdmm);
pecs = reorder(pecs);

cdmm_mean = mean(cdmm, 3);
pecs_mean = mean(pecs, 3);
cdmm_std = std(cdmm, [], 3);
pecs_std = std(pecs, [], 3);

fig{end+1} = figure("Name", "Runtime Samples");
order_label = "Order: " + string(poly_order);
type_label = ["CDMM", "PECS"];
labels = type_label + " " + order_label;
semilogy(samples, cdmm_mean, '-s'); hold on; 
semilogy(samples, pecs_mean, '-.^');
ylim([1e-4, 1e2])
ylabel("Runtime (s)");
xlabel("Number of Samples (-)")

legend(labels(:), "Location","best", "Box","off", "NumColumns", 2);

box on;

for i = 1:numel(fig)
    figure_handler.PreparePdf(fig{i}, "colorList",["#005287", "#D37611", "#00934B"]);
end
figure_handler.GenerateOutputs();
clear ans;