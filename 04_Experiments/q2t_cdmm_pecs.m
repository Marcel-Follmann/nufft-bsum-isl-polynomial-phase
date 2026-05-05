function q2t_cdmm_pecs(nr_sample_vec, poly_order_vec)
    arguments
        nr_sample_vec (:,1) double {mustBeInteger} = 32:32:256, 
        poly_order_vec (:,1) double {mustBeInteger} = 3:5
    end
    if isscalar(dbstack), close all; clc; end
    fprintf("\nQuality budget runtime anylysis: CDMM vs. PECS\n");

    max_time = 100;
    settings.nr_sample_vec = nr_sample_vec; 
    settings.poly_order_vec = poly_order_vec;
    
    data_logger = Datalogger.ChunkBased(settings,"location", "02_Data\q2t_cdmm_pecs\");
        
    for idx_samples = 1:numel(nr_sample_vec)
        for idx_poly_order = 1:numel(poly_order_vec)
            % fprintf("Samples: %d, Order: %d", nr_sample_vec(idx_samples), poly_order_vec(idx_poly_order));
            
            % Prepare Methods
            nr_samples = nr_sample_vec(idx_samples);
            poly_order = poly_order_vec(idx_poly_order);
            params = Algorithm.Params(...
                "nr_samples", nr_samples, ...
                "nr_params", poly_order, ...
                "max_mm_iters", nr_sample_vec(idx_samples)^3 ...
                );
    
            isl_pecs = 0;

            [x0, a0] = params.get_random_sequence();

            cdmm = Algorithm.CDMM.General(params);
            tic; 
            a = cdmm.run(a0);
            time_cdmm = toc;
            isl_cdmm = params.compute_isl_dB(a);
            
            pecs = Algorithm.PECS(params);
            
            a = a0;
            x = x0;
            time_pecs = 0;
            for k = 1:params.max_mm_iters
                tmp = isl_pecs;
                tic; 
                xTmp = pecs.mpslStep(x);

                [x, a] = pecs.perfectSquare(xTmp, a);
                time_pecs_tmp = toc;
                time_pecs = time_pecs + time_pecs_tmp;
                isl_pecs = params.compute_isl_dB(a);

                if k > 1
                    if abs(isl_pecs - tmp) < params.tol
                        break;
                    end
                end
                if isl_pecs <= isl_cdmm - 3 % Account for 3 dB Tolerance
                    break;
                end
                if time_pecs >= max_time
                    break;
                end
            end
            data_logger.append(nr_samples, poly_order, time_cdmm, isl_cdmm, time_pecs, isl_pecs);
            fprintf(".");
        end
    end
    data_logger.save();

    if isscalar(dbstack), keyboard, end
end


