close all; clearvars; clc; 

nr_monte_carlo = 500 - 146;
nr_samples(:,1) = round(2.^(3:0.25:6));%round(2.^(3:0.25:8));
poly_order(:,1) = 3:5;

for i = 1:nr_monte_carlo
    cd_improvement_cdmm();
    % t2q_cdmm_pecs(nr_samples, poly_order);
    % q2t_cdmm_pecs(nr_samples, poly_order);
    % rt_cdmm_pecs(nr_samples, poly_order);
    % rt_cdmm_fisl(nr_samples);
end