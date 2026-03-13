close all; clearvars; clc; 

n_samples = 2^4;
n_poly_order = 3;
n_max_iters = 10*n_samples^3;
tolerance = 1e-12;

params = Params(...
    "nr_samples", n_samples, ...
    "nr_params",n_poly_order, ...
    "max_mm_iters",n_max_iters, ...
    "tol", tolerance, ...
    "nr_candidates",1 ...
    );

rng shuffle

parfor i = 1:1e4
    [~, a0] = params.get_random_sequence();
    a(i) = a0(3);
end
histogram(a);