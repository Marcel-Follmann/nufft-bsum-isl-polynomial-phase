function mi_cdmm(nr_monte_carlo, nr_samples_vec, poly_order_vec)
    arguments
        nr_monte_carlo (1,1) double {mustBeInteger} = 150, 
        nr_samples_vec (:,1) double {mustBeInteger} = [4;16;64], 
        poly_order_vec (:,1) double {mustBeInteger} = 3:6,
    end
    if isscalar(dbstack), close all; clc; end
    fprintf("Mutual Interference Evaluation of CDMM\n");

    settings.nr_sample_vec = nr_samples_vec; 
    settings.poly_order_vec = poly_order_vec;
    
    for k = 1:numel(nr_samples_vec)
        for i = 1:numel(poly_order_vec)
            data_logger = Datalogger.ChunkBased(settings, "location", "02_Data\mi\");
            params = Algorithm.Params(...
                "nr_samples", nr_samples_vec(k), ...
                "nr_params",poly_order_vec(i) ...
                );
            
            cdmm = Algorithm.CDMM.General(params);
            
            for mc = 1:nr_monte_carlo
                [~, a0] = params.get_random_sequence();
                a(:,mc) = cdmm.run(a0);
                x(:,mc) = params.generate_signal(a(:,mc));
            end
            
            poly_order = poly_order_vec(i)-1;
            nr_samples = sprintf("N = %d",nr_samples_vec(k));
            for mc = 1:nr_monte_carlo
                for j = mc+1:nr_monte_carlo
                    y = xcorr(x(:,mc), x(:,j))/(nr_samples_vec(k));
                    % plot(20*log10(abs(y))); hold on; 
                    mi = max((abs(y)));
                    data_logger.append(mi, poly_order, nr_samples);
                end
            end
            data_logger.save();
            
            
            clear a x
            fprintf("Order: %d, Length: %d, Done.\n", poly_order_vec(i), nr_samples_vec(k));
        end
    end
    

    if isscalar(dbstack), keyboard, end
end
