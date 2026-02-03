close all; clearvars; clc; 

addpath(".\");
addpath("02_Data\");
addpath("01_Figures\");
addpath("99_Utilities\utils\");

figure_handler = FigureHandler("savePath","01_Figures\");
fig = {};

[file, path] = uigetfile('*.mat', "Select Data", "02_Data");
load("02_Data\runtime_analysis.mat");

cdmm = reorder(cdmm, nr_monte_carlos, samples, poly_order);
pecs = reorder(pecs, nr_monte_carlos, samples, poly_order);

cdmm_mean = mean(cdmm, 3);
pecs_mean = mean(pecs, 3);
cdmm_std = std(cdmm, [], 3);
pecs_std = std(pecs, [], 3);

fig{end+1} = figure("Name", "Runtime Samples");
order_label = "Order: " + string(poly_order);
type_label = ["CDMM", "PECS"];
labels = type_label + " " + order_label;
plot(samples, cdmm_mean); hold on; 
plot(samples, pecs_mean);
ylabel("Runtime (s)");
xlabel("Number of Samples (-)")

legend(labels(:), "Location","best");

box on; grid on; 

for i = 1:numel(fig)
    figure_handler.PrepareEmf(fig{i});
end
clear ans;





function x = reorder(x, nr_monte_carlos, samples, poly_order)
    x = permute( ...
        reshape(x, nr_monte_carlos,  numel(samples), numel(poly_order)), ...
        [3,2,1]);
end