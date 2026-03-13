close all; clearvars; clc; 

addpath(".\");
addpath("02_Data\");
addpath("01_Figures\");
addpath("99_Utilities\utils\");

figure_handler = FigureHandler("savePath","01_Figures\");
fig = {};

file = uigetfile('*.mat', "Select Data", "02_Data");
load("02_Data\" + string(file));

cdmm = reorder(cdmm, nr_monte_carlos, samples, poly_order);
pecs = reorder(pecs, nr_monte_carlos, samples, poly_order);

cdmm = permute(cdmm(end, end, :),[3,1,2]);
pecs = permute(pecs(end, end, :),[3,1,2]);

params = Params("nr_params",poly_order, "nr_samples",samples);

isl_cdmm = params.compute_isl_dB(permute(a_cdmm,[1,3,2])).';
isl_pecs = params.compute_isl_dB(permute(a_pecs,[1,3,2])).';

fig{end+1} = figure("Name", "Scatter Comparison");
plot(cdmm, isl_cdmm, '.', "DisplayName","CDMM"); hold on ;
plot(pecs, isl_pecs, '.', "DisplayName","PECS"); 
box on; grid on; 
xlabel("Runtime (s)");
ylabel("ISL (dB)");
legend("Location","best");

for i = 1:numel(fig)
    figure_handler.PrepareEmf(fig{i}, "colorList",["#005287", "#D37611", "#00934B"]);
end
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