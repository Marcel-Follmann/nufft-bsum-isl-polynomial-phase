close all; clearvars; clc; 




addpath(".\");
addpath("02_Data\");
addpath("01_Figures\");
addpath("99_Utilities\utils\");

figure_handler = FigureHandler("savePath","01_Figures\");
fig = {};

file = uigetfile('*.mat', "Select Data", "02_Data");
load("02_Data\" + string(file));

isl_cdmm = reorder(isl_cdmm, nr_monte_carlos, samples, poly_order);
isl_pecs = reorder(isl_pecs, nr_monte_carlos, samples, poly_order);

time_cdmm = reorder(time_cdmm, nr_monte_carlos, samples, poly_order);
time_pecs = reorder(time_pecs, nr_monte_carlos, samples, poly_order);
isl_pecs(time_pecs > 10) = NaN;
time_pecs(time_pecs > 10) = NaN;

time_pecs = mean(time_pecs,3, "omitmissing");
isl_pecs = mean(isl_pecs,3, "omitmissing");
time_cdmm = mean(time_cdmm,3, "omitmissing");
isl_cdmm = mean(isl_cdmm,3, "omitmissing");


fig{end+1} = figure("Name","Quality to Time");
semilogy(samples, time_cdmm, '-^'); hold on; 
semilogy(samples, time_pecs, '-.s');
ylim([1e-4, 1e2]);
box on; %grid on; 

labels = ["CDMM, Order: " + num2str(poly_order); "PECS, Order: " + num2str(poly_order)];
legend(labels, "Location","best", "Box","off", "NumColumns",2);
xlabel("Nr of Samples (-)");
ylabel("Runtime ($s$)");

for i = 1:numel(fig)
    figure_handler.PreparePdf(fig{i}, "colorList",["#005287", "#D37611", "#00934B"]);
end
figure_handler.GenerateOutputs();
clear ans;







function x = reorder(x, nr_monte_carlos, samples, poly_order)
    x = permute( ...
        reshape(x, nr_monte_carlos,  numel(poly_order), numel(samples)), ...
        [3,2,1]);
end