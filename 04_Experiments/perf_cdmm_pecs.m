function perf_cdmm_pecs(nr_samples, poly_order, nr_max_iters, tolerance)
    arguments
        nr_samples (1,1) double {mustBeInteger} = 128, 
        poly_order (1,1) double {mustBeInteger} = 3,
        nr_max_iters (1,1) double {mustBeInteger} = 1e5,
        tolerance (1,1) double = 1e-8
    end
    if isscalar(dbstack), close all; clc; end
    fprintf("Performance Evaluation of CDMM and PECS\n");

    rng("shuffle");
    
    params = Algorithm.Params(...
        "nr_samples", nr_samples, ...
        "nr_params",poly_order, ...
        "max_mm_iters",nr_max_iters, ...
        "tol", tolerance ...
        );
    
    pecs = Algorithm.PECS(params);
    cdmm = Algorithm.CDMM.General(params);
    
    [x0, a0] = params.get_random_sequence();
    %a0 = mod(a0, 2*pi) - pi; % as this is periodic
    
    fprintf("Run PECS\n")
    pecs.run(x0, a0, ...
        "data_logger", Datalogger.Create("file_location", "02_Data\perf_2","file_name","perf_pecs"));
    
    fprintf("Run CDMM\n")
    cdmm.run(a0, ...
        "data_logger", Datalogger.Create("file_location", "02_Data\perf_2","file_name","perf_cdmm"));

    if isscalar(dbstack), keyboard, end
end


