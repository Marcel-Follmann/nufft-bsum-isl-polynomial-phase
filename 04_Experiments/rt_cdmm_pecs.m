function rt_cdmm_pecs(nr_sample_vec, poly_order_vec)
    arguments
        nr_sample_vec (:,1) double {mustBeInteger} = 32:32:256, 
        poly_order_vec (:,1) double {mustBeInteger} = 3:5
    end
    if isscalar(dbstack), close all; clc; end
    fprintf("\nRuntime anylysis: CDMM vs. PECS\n");
    
    settings.nr_sample_vec = nr_sample_vec; 
    settings.poly_order_vec = poly_order_vec;

    data_logger = Datalogger.ChunkBased(settings, "location","02_Data\rt_cdmm_pecs\");

    for idx_samples = 1:numel(nr_sample_vec)
        for idx_poly_order = 1:numel(poly_order_vec)
            % fprintf("Samples: %d, Order: %d", nr_sample_vec(idx_samples), poly_order_vec(idx_poly_order));
    
            % Prepare Methods
            nr_samples = nr_sample_vec(idx_samples);
            poly_order = poly_order_vec(idx_poly_order);

            params = Algorithm.Params(...
                "nr_samples",nr_samples, ...
                "nr_params", poly_order ...
                );
            pecs = Algorithm.PECS(params);
            cdmm = Algorithm.CDMM.General(params);
    
            [x0, a0] = params.get_random_sequence();

            % PECS
            tic; 
            pecs.run(x0, a0);
            timePECS = toc;

            % CDMM
            tic; 
            cdmm.run(a0);
            timeCDMM = toc;
    
            data_logger.append(nr_samples, poly_order, timePECS, timeCDMM);

            fprintf(".");
        end
    end
    data_logger.save();

    if isscalar(dbstack), keyboard, end
end