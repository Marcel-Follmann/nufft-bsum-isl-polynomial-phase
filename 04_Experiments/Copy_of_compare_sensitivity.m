close all; clearvars; clc; 

addpath("99_Utilities\utils\");
fh = FigureHandler("savePath","01_Figures\");
f = {};

vec_poly_order = 3:6;
nr_samples = 4;





for i = 1:numel(vec_poly_order)
    data_logger(i) = Datalogger.Create("file_name","compare_sensitivity" + num2str(i));
    data_logger(i).log(...
    "samples", nr_samples, ...
    "poly_order", vec_poly_order(i) ...
    );

    params = Params(...
        "nr_samples",nr_samples, ...
        "nr_params", vec_poly_order(i), ...
        "tol", 1e-6, ...
        "nr_candidates", 1 ...
        );

    cdmm = CDMM.Type3(params);

    a0 = zeros(vec_poly_order(i), 1);%[x0, a0] = params.get_random_sequence();
    a_cdmm = cdmm.run(a0, "data_logger",data_logger(i), "a_range",[0, 2*pi]);
    data_logger(i).delete();
end




for i = 1:numel(vec_poly_order)
    p = vec_poly_order(i)-1;
    if isprime(p)
        mirror = linspace(0,2*pi, p+1);
    else
        mirror = linspace(0,2*pi, 3);
    end
    f{end+1} = figure("Name","Periodicity " + num2str(i));
    load("02_Data\compare_sensitivity"+num2str(i)+".mat");
    a = linspace(mirror(1),mirror(2),1e4);
    r(i,:) = 20*log10(abs(isl(end).interpolate(a)));
    plot(a, r(i,:));
    % polarplot(a, r(i,:)); hold on; 
    % rlim([floor(min(r(i,:), [],'all')/10)*10, ceil(max(r(i,:), [],'all')/10)*10]);
    % 
    % polarplot(mirror.', rlim.'.*S);
end



for i = 1:numel(f)
    fh.PrepareEmf(f{i});
end
fh.GenerateOutputs();
