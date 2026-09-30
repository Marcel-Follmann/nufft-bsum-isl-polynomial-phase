function grid_size_cdmm(nr_sample_vec, poly_order_vec, grid_size_vec)
    arguments
        nr_sample_vec (:,1) double {mustBeInteger} = 32, 
        poly_order_vec (:,1) double {mustBeInteger} = 3:4
        grid_size_vec (:,1) double = 2 * ceil(round(logspace(0,5,6)) / 2) - 1;
    end
    if isscalar(dbstack), close all; clc; end
    
    fprintf("\nGrid size influence\n");

    settings.nr_samples_vec = nr_sample_vec;
    settings.poly_order_vec = poly_order_vec;
    settings.grid_size_vec = grid_size_vec;
    data_logger = Datalogger.ChunkBased(settings, "location","02_Data\grid_size\");

    for idx_samples = 1:numel(nr_sample_vec)
        nr_samples = nr_sample_vec(idx_samples); 

        for idx_poly_order = 1:numel(poly_order_vec)
            poly_order = poly_order_vec(idx_poly_order);

            a0 = rand([1,poly_order])*2*pi;
            
            % a0 = zeros([1,poly_order])*2*pi;
            % a0(3) = rand([1,1])*2*pi;

            for idx_grid_size = 1:numel(grid_size_vec)
                grid_size = grid_size_vec(idx_grid_size);
                params = Algorithm.Params(...
                    "nr_samples", nr_samples, ...
                    "nr_params", poly_order, ...
                    "nr_candidates", grid_size, ...
                    "max_mm_iters", nr_sample_vec(idx_samples)^3, ...
                    "nr_cd_iterations", 1000 ...
                    );
        
                % [~, a0] = params.get_random_sequence();
                cdmm = Algorithm.CDMM.General(params);
                tic; 
                a = cdmm.run(a0);
                time = toc; 
                x = params.generate_signal(a);
                isl = params.compute_isl_dB(x);
    
                data_logger.append(nr_samples, poly_order, grid_size, isl, time);
                fprintf(".");
            end
        end
    end
    data_logger.save();

    if isscalar(dbstack), keyboard, end
end