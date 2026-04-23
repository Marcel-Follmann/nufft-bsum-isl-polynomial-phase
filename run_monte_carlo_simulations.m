close all; clearvars; clc; 

nr_monte_carlo = 1000;
nr_samples(:,1) = round(2.^(3:0.25:8));%8:8:32:256;
poly_order(:,1) = 3:4:5;

rt_cdmm_pecs(nr_monte_carlo, nr_samples, poly_order);
t2q_cdmm_pecs(nr_monte_carlo, nr_samples, poly_order);
q2t_cdmm_pecs(nr_monte_carlo, nr_samples, poly_order);
rt_cdmm_fisl(nr_monte_carlo, nr_samples);