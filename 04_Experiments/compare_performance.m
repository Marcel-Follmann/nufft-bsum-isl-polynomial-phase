close all; clearvars; clc; 

addpath(".\");
addpath("02_Data\");
addpath("99_Utilities\utils\");
rng("shuffle");


n_samples = 2^4;
n_poly_order = 6;
n_max_iters = n_samples^3;
tolerance = 1e-8;

pecs = PECS("nr_samples", n_samples, "nr_params",n_poly_order, "save_data", true, "file_name", "pecs_performance", ...
    "max_mm_iters",n_max_iters, "tol", tolerance);
cdmm = CDMMT3("nr_samples", n_samples, "nr_params",n_poly_order, "save_data", true, "file_name", "cdmm_performance", ...
    "max_mm_iters",n_max_iters, "tol", tolerance, "nr_candidates",1);

[x0, a0] = pecs.getRandomSequence();
a0 = mod(a0, 2*pi); % as this is periodic
pecs.mpsl(x0, a0);
cdmm.type2(a0);

delete(cdmm);
delete(pecs);


