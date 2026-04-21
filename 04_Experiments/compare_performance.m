close all; clearvars; clc; 

rng("shuffle");

n_samples = 100;
n_poly_order = 3;
n_max_iters = 1e5;
tolerance = 1e-7;

params = Algorithm.Params(...
    "nr_samples", n_samples, ...
    "nr_params",n_poly_order, ...
    "max_mm_iters",n_max_iters, ...
    "tol", tolerance, ...
    "nr_candidates",501 ...
    );

pecs = Algorithm.PECS(params);
cdmm = Algorithm.CDMM.General(params);

[x0, a0] = params.get_random_sequence();
a0 = mod(a0, 2*pi); % as this is periodic

pecs.run(x0, a0, ...
    "data_logger", Datalogger.Create("file_name","pecs_performance"));

cdmm.run(a0, ...
    "data_logger", Datalogger.Create("file_name","cdmm_performance"));


