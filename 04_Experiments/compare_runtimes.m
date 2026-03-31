close all; clearvars; clc; 

addpath(".\");
addpath("02_Data\");
addpath("99_Utilities\utils\");


nr_monte_carlo = 1;
nr_samples(:,1) =  150;%round(2.^(4:.125:8));%2.^(4:1:5); % max 8
poly_order(:,1) = 3;%(3:5); % max 4
tol = 1e-4;
relErrorTol = 1e-2; 

data_logger = Datalogger.Create("file_name","runtime_analysis", "add_timestamp",true);
data_logger.log(...
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
            "max_mm_iters", nr_samples(idx_samples)^3, ...
            "tol", tol, ...
            "nr_candidates",1 ...
            );
        pecs = PECS(params);
        cdmm = CDMM.Type3(params);

        relError = Inf; 
        timePECS = [];
        timeCDMM = [];
        i = 0;
        while any([relError > relErrorTol, i < 10])
            i = i+1;
            [x0, a0] = params.get_random_sequence();

            % PECS
            tic; 
            a_pecs= pecs.run(x0, a0);
            timePECS(i) = toc;

            % CDMM
            tic; 
            a_cdmm = cdmm.run(a0);
            timeCDMM(i) = toc;

            t = tinv(0.95, i);
            sigma(1) = sqrt(var(timeCDMM))/i;
            sigma(2) = sqrt(var(timePECS))/i;
            x(1) = mean(timeCDMM);
            x(2) = mean(timePECS);
            relError = 2*sigma*t./x;

            fprintf(".");
        end

        data_logger.log("pecs", x(2), "cdmm", x(1), "nr_monte_carlos", i);
        fprintf("\n");
    end
end
data_logger.delete();