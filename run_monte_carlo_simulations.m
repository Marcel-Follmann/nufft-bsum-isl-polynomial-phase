close all; clearvars; clc; 

nr_monte_carlo = 5;
nr_samples(:,1) = 32:34;%:32:256;
poly_order(:,1) = 3;%:5;

rt_cdmm_pecs(nr_monte_carlo, nr_samples, poly_order);
t2q_cdmm_pecs(nr_monte_carlo, nr_samples, poly_order);
q2t_cdmm_pecs(nr_monte_carlo, nr_samples, poly_order);
rt_cdmm_fisl(nr_monte_carlo, nr_samples);