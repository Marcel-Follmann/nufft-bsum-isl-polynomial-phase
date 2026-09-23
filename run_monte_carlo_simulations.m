close all; clearvars; clc; 

nr_monte_carlo = 500;
nr_samples(:,1) = round(2.^(3:0.25:6));%round(2.^(3:0.25:8));
poly_order(:,1) = 3:5;


for i = 1:nr_monte_carlo
    
    % isl_bound(nr_samples, poly_order);
    % cd_improvement_cdmm();
    t2q_cdmm_pecs(nr_samples, poly_order);
    % q2t_cdmm_pecs(nr_samples, poly_order);
    % rt_cdmm_pecs(nr_samples, poly_order);
    % rt_cdmm_fisl(nr_samples);
end
% clearvars; 
% nr_monte_carlo = 500;
% nr_samples(:,1) = [32];%round(2.^(3:0.25:8));
% poly_order(:,1) = [3, 4];
% grid_size = 2 * ceil(round(logspace(0,5,6)) / 2) - 1;
% for i = 1:nr_monte_carlo
%     grid_size_cdmm(nr_samples, poly_order, grid_size)
% end
% 
