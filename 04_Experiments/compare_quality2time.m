close all; clearvars; clc; 

nr_monte_carlo = 50;
nr_samples(:,1) = 16:8:64;%2.^(4:7); % max 8
poly_order(:,1) = (3:5); % max 4
tol = 1e-4;
max_time = Inf;

data_logger = Datalogger.Create("file_name","quality2time_analysis", "add_timestamp",true);
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
            "max_mm_iters", nr_samples(idx_samples)^3, ...
            "tol", tol, ...
            "nr_candidates",1 ...
            );

        isl_pecs = zeros(nr_monte_carlo,1);
        for i = 1:nr_monte_carlo
            [x0, a0] = params.get_random_sequence();

            cdmm = CDMM.Type3(params);
            tic; 
            a = cdmm.run(a0);
            time_cdmm(i) = toc;
            isl_cdmm(i) = params.compute_isl_dB(a);
            
            pecs = PECS(params);
            
            a = a0;
            x = x0;
            time_pecs(i) = 0;
            for k = 1:1e6
                tmp = isl_pecs(i);
                tic; 
                xTmp = pecs.mpslStep(x);

                [x, a] = pecs.perfectSquare(xTmp, a);
                time_pecs_tmp = toc;
                time_pecs(i) = time_pecs(i) + time_pecs_tmp;
                isl_pecs(i) = params.compute_isl_dB(a);

                if i > 1
                    if abs(isl_pecs(i) - tmp) < 1e-6
                        break;
                    end
                end
                if isl_pecs(i) < isl_cdmm
                    break;
                end
                if time_pecs(i) >= min([100*time_cdmm(i),max_time])
                    break;
                end
            end
            fprintf(".");
        end
        data_logger.log("time_cdmm", time_cdmm, "isl_cdmm", isl_cdmm, "time_pecs", time_pecs, "isl_pecs", isl_pecs);
        fprintf("\n");
    end
end

clearvars;


