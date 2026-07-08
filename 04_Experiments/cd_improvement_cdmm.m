function cd_improvement_cdmm(nr_samples_vec, poly_order)
    arguments
        nr_samples_vec (:,1) double {mustBeInteger} = [16, 32, 64], 
        poly_order (1,1) double {mustBeInteger} = 5
    end
    if isscalar(dbstack), close all; clc; end
    fprintf("Coordinate Improvement of CDMM\n");

    rng("shuffle");
    
    % fh = Figurehandler.Create("savePath","01_Figures\");
    % f = {};
    % colorlist = flip(["#00934B", "#D37611", "#005287"]);
    settings.nr_sample_vec = nr_samples_vec; 
    settings.poly_order = poly_order;
    data_logger_external = Datalogger.ChunkBased(settings, "location", "02_Data\ci\");

    for k = 1:numel(nr_samples_vec)
        params = Algorithm.Params(...
        "nr_samples", nr_samples_vec(k), ...
        "nr_params",poly_order ...
        );
        
        cdmm = Algorithm.CDMM.General(params);
        data_logger_internal = Datalogger.Create("file_location", "02_Data\ci", "file_name","ci");
        [~, a0] = params.get_random_sequence();
        cdmm.run(a0, "data_logger", data_logger_internal);
        data_logger_internal.delete();
        while ~exist("02_Data\ci\ci.mat", "file"), pause(0.001); end
        load("02_Data\ci\ci.mat");
        delete("02_Data\ci\ci.mat");

        for i = 1:size(a,2)
            x = params.generate_signal(a(:,i));
            isl_coordinate(i) = params.compute_isl_dB(x);
            % data_logger_external.append(a, coordinate, isl);
        end
        
        
        tmp = diff(isl_coordinate);
        for i = 1:size(a,2)-1
            coordinate = i + 1;
            nr_samples = "N = " + num2str(nr_samples_vec(k));
            improvement = tmp(i);
            data_logger_external.append(improvement, coordinate, nr_samples);
        end
    end
    data_logger_external.save();
    fprintf(".");

    
    
    % f{end+1} = figure("Name", "ISL Improvement");
    % group = categorical(T.samples);
    % group = reordercats(group, {'N = 16','N = 32', 'N = 64'});
    % b = boxchart(T.coordinate, -T.improvement, "GroupByColor",T.samples, "Notch","off");
    % ylabel("ISL Improvement (dB)")
    % xlabel("Coordinate (-)");
    % 
    % ax = gca;
    % ds = 1/2;
    % xticks(2:poly_order-1)
    % xlim([2 - ds, poly_order + ds-1]);
    % 
    % for i = 1:numel(b)
    %     b(i).JitterOutliers = 'on';
    %     b(i).MarkerStyle = '.';
    %     b(i).BoxFaceColor = hex2rgb(colorlist(i));
    %     b(i).MarkerColor = hex2rgb(colorlist(i));
    % end
    % lgd = legend("Location","best", "Box", "off");
    % 
    % fh.PreparePdf(f{1});
    % ax.XMinorTick = 'off';
    % ax.YGrid ="on";
    % fh.GenerateOutputs();

    if isscalar(dbstack), keyboard, end
end



