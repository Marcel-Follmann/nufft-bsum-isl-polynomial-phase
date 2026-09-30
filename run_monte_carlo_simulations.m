close all; clearvars; clc; 

%%
nr_monte_carlo = 99;
for i = 1:nr_monte_carlo
    bit_res_cdmm();
    % grid_size_cdmm();
end

%%
nr_monte_carlo = 0; 
for i = 1:nr_monte_carlo
    cd_improvement_cdmm();
end

%%
nr_monte_carlo = 0;
nr_samples(:,1) = round(2.^(3:0.25:8)); 
poly_order(:,1) = 3:5; 

timer = tic; 
for i = 1:nr_monte_carlo
    rt_cdmm_pecs(nr_samples, poly_order);
    rt_cdmm_fisl(nr_samples);
    t2q_cdmm_pecs(nr_samples, poly_order);
end
toc(timer);
