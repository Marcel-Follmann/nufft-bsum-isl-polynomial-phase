function sen_cdmm(nr_samples, poly_order)
    arguments
        nr_samples (:,1) double {mustBeInteger} = 2:128,%round(linspace(2,100,50)), 
        poly_order (:,1) double {mustBeInteger} = 3,
    end
    if isscalar(dbstack), close all; clc; end
    fprintf("Sensitivity Evaluation of CDMM\n");
    
    data_logger = Datalogger.Create("file_name","sen");
    
    for i = 1:numel(nr_samples)
        data_logger.log(...
        "samples", nr_samples(i), ...
        "poly_order", poly_order ...
        );
    
        params = Algorithm.Params(...
            "nr_samples",nr_samples(i), ...
            "nr_params", poly_order ...
            );
    
        cdmm = Algorithm.CDMM.General(params);
    
        a0 = zeros(poly_order, 1);%[x0, a0] = params.get_random_sequence();
        a_cdmm = cdmm.run(a0, "data_logger",data_logger, "a_range",[0, 2*pi]);
        fprintf(".");
    end
    fprintf("\n");
    data_logger.delete();

    if isscalar(dbstack), keyboard, end
end


