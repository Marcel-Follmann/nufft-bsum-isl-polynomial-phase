function rt_cdmm_pecs(nr_monte_carlo, nr_samples, poly_order)
    arguments
        nr_monte_carlo (1,1) double {mustBeInteger} = 10, 
        nr_samples (:,1) double {mustBeInteger} = 32:32:256, 
        poly_order (:,1) double {mustBeInteger} = 3:5
    end
    if isscalar(dbstack), close all; clc; end
    fprintf("Runtime anylysis: CDMM vs. PECS\n");
    
    data_logger = Datalogger.Create( ...
        "file_location", "02_Data\rt_cdmm_pecs", ...
        "file_name", "MC_"+num2str(nr_monte_carlo)+"_rt_cdmm_pecs", ...
        "add_timestamp", true ...
        );
    
    data_logger.log(...
        "samples", nr_samples, ...
        "poly_order", poly_order, ...
        "nr_monte_carlos", nr_monte_carlo ...
        );
    
    for idx_samples = 1:numel(nr_samples)
        for idx_poly_order = 1:numel(poly_order)
            fprintf("Samples: %d, Order: %d -> ", nr_samples(idx_samples), poly_order(idx_poly_order));
    
            % Prepare Methods
            params = Algorithm.Params(...
                "nr_samples",nr_samples(idx_samples), ...
                "nr_params", poly_order(idx_poly_order) ...
                );
            pecs = Algorithm.PECS(params);
            cdmm = Algorithm.CDMM.General(params);
    
            for i = 1:nr_monte_carlo
                [x0, a0] = params.get_random_sequence();
    
                % PECS
                tic; 
                a_pecs= pecs.run(x0, a0);
                timePECS = toc;
    
                % CDMM
                tic; 
                a_cdmm = cdmm.run(a0);
                timeCDMM = toc;
    
                data_logger.log("pecs", timePECS, "cdmm", timeCDMM);
                if mod(i, nr_monte_carlo/100) == 0
                    fprintf(".");
                end
            end
            fprintf("\n");
        end
    end

    if isscalar(dbstack), keyboard, end
end