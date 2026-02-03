close all; clearvars; clc; 

addpath(".\");
addpath("02_Data\");

data_logger = DataLogger("log_enabled",true, "file_name","runtime_analysis");

nMonteCarlo = 50;
nr_samples(:,1) = 10:10:100;
poly_order(:,1) = 3:4;

data_logger.log(...
    "nr_monte_carlos", nMonteCarlo, ...
    "samples", nr_samples, ...
    "poly_order", poly_order ...
    );

for idx_samples = 1:numel(nr_samples)
    for idx_poly_order = 1:numel(poly_order)
        fprintf("Samples: %d, Order: %d -> ", nr_samples(idx_samples), poly_order(idx_poly_order));

        % Prepare Methods
        pecs = PECS("nr_samples",nr_samples(idx_samples), "nr_params",poly_order(idx_poly_order));
        cdmm = CDMM("nr_samples",nr_samples(idx_samples), "nr_params", poly_order(idx_poly_order), "nr_candidates",1);

        for iMonteCarlo = 1:nMonteCarlo
            [x0, a0] = pecs.getRandomSequence();

            % PECS
            tic; 
            pecs.mpsl(x0, a0);
            
            timePECS = toc;

            % CDMM
            tic; 
            cdmm.type2(a0);
            timeCDMM = toc;

            data_logger.log("pecs", timePECS, "cdmm", timeCDMM);

            fprintf(".");
        end
        fprintf("\n");
    end
end

