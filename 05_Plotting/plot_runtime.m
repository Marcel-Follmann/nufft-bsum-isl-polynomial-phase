close all; clearvars; clc; 

figure_handler = Figurehandler.Create("savePath","01_Figures\");
fig = {};

file = uigetfile('*.mat', "Select Data", "02_Data");
load("02_Data\" + string(file));

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
semilogy(samples, cdmm_mean, '-s'); hold on; 
semilogy(samples, pecs_mean, '-.^');
ylim([1e-4, 1e2])
ylabel("Runtime (s)");
xlabel("Number of Samples (-)")
% x = [samples; flip(samples)];
% y = [cdmm_mean - 3*cdmm_std; flip(cdmm_mean + 3*cdmm_std,1)];
% fill(x,y(:,1), hex2rgb("#005287"), "FaceAlpha", 0.1, "EdgeAlpha", 0)
% fill(x,y(:,2), hex2rgb("#00934B"), "FaceAlpha", 0.1, "EdgeAlpha", 0)
% y = [pecs_mean - 3*pecs_std; flip(pecs_mean + 3*pecs_std,1)];
% fill(x,y(:,1), hex2rgb("#D37611"), "FaceAlpha", 0.1, "EdgeAlpha", 0)
% fill(x,y(:,2), hex2rgb("#005287"), "FaceAlpha", 0.1, "EdgeAlpha", 0)

legend(labels(:), "Location","best", "Box","off", "NumColumns", 2);

box on; %grid on; 

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

function rgb = hex2rgb(hex)
	% Remove the # symbol if it's there
    hex = strrep(hex, '#', '');
	hex=char(hex);
	
	% Convert hex values to decimal and normalize by dividing by 255
    r = hex2dec(hex(1:2)) / 255;
    g = hex2dec(hex(3:4)) / 255;
    b = hex2dec(hex(5:6)) / 255;

    % Output RGB triplet
    rgb = [r, g, b];
end