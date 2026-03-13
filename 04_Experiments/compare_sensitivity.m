close all; clearvars; clc; 

addpath("99_Utilities\utils\");
fh = FigureHandler("savePath","01_Figures\");
f = {};

poly_order = 3;
nr_samples = round(linspace(5,50,49));

data_logger = Datalogger.Create("file_name","compare_sensitivity");



for i = 1:numel(nr_samples)
    data_logger.log(...
    "samples", nr_samples(i), ...
    "poly_order", poly_order ...
    );

    params = Params(...
        "nr_samples",nr_samples(i), ...
        "nr_params", poly_order, ...
        "tol", 1e-6, ...
        "nr_candidates", 1 ...
        );

    cdmm = CDMM.Type3(params);

    a0 = zeros(poly_order, 1);%[x0, a0] = params.get_random_sequence();
    a_cdmm = cdmm.run(a0, "data_logger",data_logger, "a_range",[0, 2*pi]);
end
data_logger.delete();





load("02_Data\compare_sensitivity.mat");

a = linspace(0,pi,1e4);

for i = 1:numel(isl)
    r(i,:) = 20*log10(abs(isl(i).interpolate(a)));
end

[A,Z] = meshgrid(2*a,samples);

X = r .* cos(A);
Y = r .* sin(A);

f{end+1} = figure("Name", "3D Sensitivity");
surf(X,Y,Z)
shading interp
daspect([1 1 max(samples)/max(X, [], "all")]); %axis equal
set(gca,'ZDir','reverse')
xticks([])
yticks([])
box on;
grid on; grid minor;
view(0,90);
colorbar;
%grid off; %grid on; grid minor; 
zlim([2,max(nr_samples)])
zlabel("Sequence Length N (-)");
% xlabel("$ISL_{dB}(a, N) \cdot \cos(a)$");
% ylabel("$ISL_{dB}(a, N) \cdot \sin(a)$");

for i = 1:numel(f)
    fh.PrepareEmf(f{i});
end
fh.GenerateOutputs();
