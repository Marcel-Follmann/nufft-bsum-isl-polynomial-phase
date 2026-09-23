close all; clearvars; clc;

data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\isl_bound\");
load("02_Data\isl_bound\settings.mat");

[G_samples, samples] = findgroups(data.nr_samples);
[G_order, order] = findgroups(data.poly_order);

for i = 1:height(data)
    nr_samples = data(i, :).nr_samples;
    poly_order = data(i,:).poly_order;
    a = cell2mat(data(i,:).a);

    params = Algorithm.Params(...
        "nr_samples", nr_samples, ...
        "nr_params", poly_order, ...
        "max_mm_iters", nr_samples^3 ...
        );

    x = params.generate_signal(a);

    isl(i, G_samples(i), G_order(i)) = params.compute_isl_dB(x);
end

boxchart(isl(:,:,1)); hold on; 
boxchart(isl(:,:,2));
boxchart(isl(:,:,3));
legend("Location","best")
hold on;
