close all; clearvars; clc; 

addpath(".\");
addpath("02_Data\");
addpath("99_Utilities\utils\");


nr_monte_carlo = 5;
nr_samples(:,1) = 2.^(3:5); % max 8
poly_order(:,1) = 3:5; % max 4

n_max_iter = 1e4; 
tol = 1e-9;

data_logger = Datalogger.Create("file_name","runtime_analysis", "add_timestamp",true);
data_logger.log(...
    "nr_monte_carlos", nr_monte_carlo, ...
    "samples", nr_samples, ...
    "poly_order", poly_order ...
    );

for idx_samples = 1:numel(nr_samples)
    for idx_poly_order = 1:numel(poly_order)
        fprintf("Samples: %d, Order: %d -> ", nr_samples(idx_samples), poly_order(idx_poly_order));

        % Prepare Methods
        params = Params(...
            "nr_samples",nr_samples(idx_samples), ...
            "nr_params", poly_order(idx_poly_order), ...
            "max_mm_iters", n_max_iter, ...
            "tol", tol, ...
            "nr_candidates",1 ...
            );
        pecs = PECS(params);
        cdmm = CDMM.Type3(params);

        timeCDMM = zeros(nr_monte_carlo, 1);
        timePECS = zeros(nr_monte_carlo, 1);
        for i = 1:nr_monte_carlo
            [x0, a0] = params.get_random_sequence();

            % PECS
            tic; 
            pecs.mpsl(x0, a0);
            timePECS(i) = toc;

            % CDMM
            tic; 
            cdmm.run(a0);
            timeCDMM(i) = toc;

            fprintf(".");
        end


        data_logger.log("pecs", timePECS, "cdmm", timeCDMM);
        fprintf("\n");
    end
end
data_logger.delete();