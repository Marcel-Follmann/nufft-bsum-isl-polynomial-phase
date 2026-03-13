close all; clearvars; clc; 

addpath(".\");
addpath("02_Data\");
addpath("99_Utilities\utils\");
rng("shuffle");

n_samples = 2^4;%2^7;
n_poly_order = 3;%5;
n_max_iters = 1e2;%10*n_samples^3;
tolerance = 1e-4;

params = Params(...
    "nr_samples", n_samples, ...
    "nr_params",n_poly_order, ...
    "max_mm_iters",n_max_iters, ...
    "tol", tolerance, ...
    "nr_candidates",1 ...
    );

pecs = PECS(params);
cdmm = CDMM.Type3(params);

[x0, a0] = params.get_random_sequence();
a0 = mod(a0, 2*pi); % as this is periodic

pecs.run(x0, a0, ...
    "data_logger", Datalogger.Create("file_name","pecs_performance"));

cdmm.run(a0, ...
    "data_logger", Datalogger.Create("file_name","cdmm_performance"));


