close all; clearvars; clc; 

figure_handler = Figurehandler.Create("savePath","01_Figures\");
fig = {};


[file, location] = uigetfile('*.mat', "Select Data", string(pwd) + "\02_Data\q2t_cdmm_pecs");
load(string(location) + string(file));


reorder = @(data) reorder_data(data, [3, 2, 1], ...
    nr_monte_carlos, poly_order, samples); % Resembles for loop order in experiments file

isl_cdmm = reorder(isl_cdmm);
isl_pecs = reorder(isl_pecs);
time_cdmm = reorder(time_cdmm);
time_pecs = reorder(time_pecs);

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

