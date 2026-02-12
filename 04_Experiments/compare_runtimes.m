close all; clearvars; clc; 

addpath(".\");
addpath("02_Data\");

data_logger = DataLogger("log_enabled",true, "file_name","runtime_analysis", "add_timestamp",true);

nr_monte_carlo = 20;
nr_samples(:,1) = 2.^(3:6); % max 8
poly_order(:,1) = 3:4; % max 4

data_logger.log(...
    "nr_monte_carlos", nr_monte_carlo, ...
    "samples", nr_samples, ...
    "poly_order", poly_order ...
    );

for idx_samples = 1:numel(nr_samples)
    for idx_poly_order = 1:numel(poly_order)
        fprintf("Samples: %d, Order: %d -> ", nr_samples(idx_samples), poly_order(idx_poly_order));

        % Prepare Methods
        pecs = PECS("nr_samples",nr_samples(idx_samples), "nr_params",poly_order(idx_poly_order), "max_mm_iters", 1e4);
        cdmm = CDMMT3("nr_samples",nr_samples(idx_samples), "nr_params", poly_order(idx_poly_order), "max_mm_iters", 1e4, "nr_candidates",1);

        timeCDMM = zeros(nr_monte_carlo, 1);
        timePECS = zeros(nr_monte_carlo, 1);
        for i = 1:nr_monte_carlo
            [x0, a0] = pecs.getRandomSequence();

            % PECS
            tic; 
            pecs.mpsl(x0, a0);
            timePECS(i) = toc;

            % CDMM
            tic; 
            cdmm.type2(a0);
            timeCDMM(i) = toc;

            fprintf(".");
        end


        data_logger.log("pecs", timePECS, "cdmm", timeCDMM);
        fprintf("\n");
    end
end

data_logger.delete();