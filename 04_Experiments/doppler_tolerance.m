function doppler_tolerance(nr_samples, nr_params, nr_monte_carlo)
    arguments
        nr_samples (:,1) double {mustBeInteger} = 64, 
        nr_params (:,1) double {mustBeInteger} = 3:5,
        nr_monte_carlo (:,1) double {mustBeInteger} = 400
    end
    if isscalar(dbstack), close all; clc; end

    fprintf("\nDoppler Tolerance Analysis\n");

    rng("shuffle");
    for i = 1:numel(nr_params)
        params(i) = Algorithm.Params(...
            "nr_samples", nr_samples, ...
            "nr_params",nr_params(i) ...
        );
        cdmm(i) = Algorithm.CDMM.General(params(i));
    end
    fisl = Algorithm.FISL(params(1));

    settings.params = params;
    data_logger = Datalogger.ChunkBased(settings, "location","02_Data\dt\");

    for mc = 1:nr_monte_carlo
        for i = 1:numel(nr_params)
            [x0, a0{i}] = params(i).get_random_sequence();

            a = cdmm(i).run(a0{i});
            x_cdmm(:,i) = params(i).generate_signal(a);
            Q = pow2db(params(i).compute_doppler_tolerance(x_cdmm(:,i)));
            type = i;
            data_logger.append(type, Q);
        end
        x_fisl = fisl.run(x0);
        Q = pow2db(params(1).compute_doppler_tolerance(x_fisl));
        type = i+1;
        data_logger.append(type, Q);

        fprintf(".");
    end
    data_logger.save();
   
    if isscalar(dbstack), keyboard, end
end



