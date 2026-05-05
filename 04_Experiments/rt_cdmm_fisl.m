function rt_cdmm_fisl(nr_sample_vec)
    arguments
        nr_sample_vec (:,1) double {mustBeInteger} = 32:32:256, 
    end
    if isscalar(dbstack), close all; clc; end
    fprintf("\nRuntime anylysis: CDMM vs. FISL\n");

    settings.nr_sample_vec = nr_sample_vec;
    data_logger = Datalogger.ChunkBased(settings, "location","02_Data\rt_cdmm_fisl\");
    %% Simulation
    for i = 1:numel(nr_sample_vec)
        nr_samples = nr_sample_vec(i);
        params = Algorithm.Params(...
        "nr_samples", nr_sample_vec(i), ...
        "nr_params", 3 ...
        );
    
        % Initialization
        cdmm = Algorithm.CDMM.General(params);
        cdmm_ood = Algorithm.CDMM.OOD(params);
        fisl = Algorithm.FISL(params);
        
        % fprintf("P: %d ...", nr_sample_vec(i));
        
        [x0, a0] = params.get_random_sequence();
    
     
        tic; 
        cdmm_ood.run(a0, "a_range",a0(3)+[-.5, .5]); 
        time_cdmm_ood = toc; 

        tic; 
        cdmm.run(a0, "a_range",a0(3)+[-.5, .5]); 
        time_cdmm = toc; 

        tic; 
        fisl.run(x0); 
        time_fisl = toc; 

        data_logger.append(nr_samples, time_cdmm_ood, time_cdmm, time_fisl);
    
        fprintf(".");
    end
    data_logger.save();
    close all; 
    
    % %% Visualization
    % f{end+1} = figure("Name","Runtime FISL vs CDMM A");
    % b = boxchart(T.samples, T.time, 'GroupByColor',T.type, 'Notch','on');
    % xlabel("Sequence Length (-)");
    % ylabel("Runtime (s)");
    % ds = (nr_sample_vec(2) - nr_sample_vec(1))/2;
    % xticks(nr_sample_vec)
    % xlim([nr_sample_vec(1) - ds, nr_sample_vec(end) + ds]);
    % ax = gca;
    % ax.XMinorTick = 'off';
    % lgd = legend("Location","northwest", "Box","off");
    % lgd.IconColumnWidth = 10;
    % 
    % for i = 1:numel(b)
    %     b(i).JitterOutliers = 'on';
    %     b(i).MarkerStyle = '.';
    %     b(i).BoxFaceColor = hex2rgb(colorlist(i));
    %     b(i).MarkerColor = hex2rgb(colorlist(i));
    % end
    % 
    % fh.PreparePdf(f{1});
    % ax.YScale ="log";
    % ax.XMinorTick = 'off';
    % fh.GenerateOutputs();

   if isscalar(dbstack), keyboard, end
end



