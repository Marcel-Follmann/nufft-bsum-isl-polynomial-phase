close all; clearvars; clc; 

fh = Figurehandler.Create("savePath","01_Figures\");
f = {};

poly_order = 3;
nr_samples = round(linspace(2,100,10));

% data_logger = Datalogger.Create("file_name","compare_sensitivity");
% 
% 
% 
% for i = 1:numel(nr_samples)
%     data_logger.log(...
%     "samples", nr_samples(i), ...
%     "poly_order", poly_order ...
%     );
% 
%     params = Algorithm.Params(...
%         "nr_samples",nr_samples(i), ...
%         "nr_params", poly_order, ...
%         "tol", 1e-6, ...
%         "nr_candidates", 1 ...
%         );
% 
%     cdmm = Algorithm.CDMM.General(params);
% 
%     a0 = zeros(poly_order, 1);%[x0, a0] = params.get_random_sequence();
%     a_cdmm = cdmm.run(a0, "data_logger",data_logger, "a_range",[0, 2*pi]);
% end
% data_logger.delete();

load("02_Data\compare_sensitivity.mat");

a = linspace(-pi/2,pi/2,1e3);

for i = 1:numel(isl)
    r(i,:) = 20*log10(abs(isl(i).interpolate(a)));
end

f{end+1} = figure("Name", "Sensitivity");
% imagesc(a/pi, samples, r)
% shading interp
imagesc(a/pi, samples, r);
axis xy;
% view(0,90)
map = CreateColormap((["#005287", "#FFFFFF", "#D37611"]));
colormap(map)
c = colorbar;
xlabel("$a_2 / \pi$ ($-$)");
ylabel("Sequence Length ($-$)");
ylabel(c, "ISL (dB)");
ticks = yticks;
ticks = [2, ticks(ticks > 2)];
yticks(ticks);
ylim([ticks(1), ticks(end)]);

for i = 1:numel(f)
    fh.PreparePdf(f{i});
end
fh.GenerateOutputs();


function hex = rgb2hex(rgb)
	r = round(rgb(1)*255);
	g = round(rgb(2)*255);
	b = round(rgb(3)*255);

	hex = sprintf("#%02s%02s%02s",dec2hex(r),dec2hex(g),dec2hex(b));
end

function map = CreateColormap(colorlist)
    N = 256;

    for iColor = 1:numel(colorlist)
        rgbList(iColor, :) = hex2rgb(colorlist(iColor));
    end

    map = interp1(linspace(0,1,iColor), rgbList, linspace(0,1,N));
end
