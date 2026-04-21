close all; clearvars; clc; 

addpath("99_Utilities");
fh = Figurehandler.Create("savePath","01_Figures\");
f = {};

vec_poly_order = 3:5;
nr_samples = 6;


for i = 1:numel(vec_poly_order)
    data_logger(i) = Datalogger.Create("file_name","compare_sensitivity" + num2str(i));
    data_logger(i).log(...
    "samples", nr_samples, ...
    "poly_order", vec_poly_order(i) ...
    );

    params = Algorithm.Params(...
        "nr_samples",nr_samples, ...
        "nr_params", vec_poly_order(i), ...
        "tol", 1e-6, ...
        "nr_candidates", 1 ...
        );

    cdmm = Algorithm.CDMM.General(params);

    a0 = zeros(vec_poly_order(i), 1);%[x0, a0] = params.get_random_sequence();
    a_cdmm = cdmm.run(a0, "data_logger",data_logger(i), "a_range",[0, 2*pi]);
    data_logger(i).delete();
end




for i = 1:numel(vec_poly_order)
    p = vec_poly_order(i)-1;
    if isscalar(unique(factor(p)))
        mirror = linspace(0,2*pi, p+1);
    else
        mirror = linspace(0,2*pi, 2);
    end
    f{end+1} = figure("Name","Periodicity " + num2str(p));
    load("02_Data\compare_sensitivity"+num2str(i)+".mat");
    a1 = linspace(-mirror(2)/2,mirror(2)/2,1e4);
    a2 = linspace(mirror(2)/2,2*pi - mirror(2)/2,1e4);
    r1(i,:) = 20*log10(abs(isl(end).interpolate(a1)));
    r2(i,:) = 20*log10(abs(isl(end).interpolate(a2)));

    % plot(a/pi, r(i,:));
    polarplot(a1, r1(i,:), "DisplayName","$\Phi_\text{ISL}$"); hold on; 
    polarplot(a2, r2(i,:), "DisplayName","$\tilde{\Phi}_\text{ISL}$"); hold on; 
    rlim([floor(min(r1(i,:), [],'all')/10)*10, ceil(max(r1(i,:), [],'all')/10)*10]);
    % thetatickformat("%g")
    thetaticks([0, 90, 180, 270])
    thetaticklabels({"$0$", "$\frac{1}{2}\pi$", "$\pi$", "$\frac{3}{2}\pi$"})
    % lgd = legend("Location","best", "Box","off");
    % lgd.ItemTokenSize = [10 10];
    % polarplot(mirror.', rlim.'.*S);
    fprintf("Nr. Samples: %d, Poly. Order: %d, Min. ISL: %.2f dB, Max. ISL: %.2f dB\n", nr_samples, vec_poly_order(i), min(r1(i,:)), max(r1(i,:)));
end

width = 12;%4.4;

for i = 1:numel(f)
    fh.PrepareSvg(f{i}, "width", width, "height",width);
end
fh.GenerateOutputs();
