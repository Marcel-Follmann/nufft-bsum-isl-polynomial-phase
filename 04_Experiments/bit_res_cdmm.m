function bit_res_cdmm(nr_samples, nr_params)
    arguments
        nr_samples (:,1) double {mustBeInteger} = 64, 
        nr_params (:,1) double {mustBeInteger} = 4
    end
    if isscalar(dbstack), close all; clc; end
    fprintf("\nBit Resolution anylysis\n");
    
    params = Algorithm.Params(...
        "nr_samples", nr_samples, ...
        "nr_params", nr_params, ...
        "tol", 5e-6, ...
        "nufft_e", 1e-6 ...
    );
    settings.params = params; 
    data_logger = Datalogger.ChunkBased(settings, "location","02_Data\bit_res_3\");
    %% Simulation
    % Initialization
    cdmm = Algorithm.CDMM.General(params);
    
    % a0 = 2*pi*rand(nr_params,1);
    [~, a0] = params.get_random_sequence();

    a = cdmm.run(a0).'; 
    data_logger.append(a);

    data_logger.save();
    close all; 
   
   if isscalar(dbstack), keyboard, end
end



