function plt_rt_cdmm_pecs(fh)
    arguments
        fh (1,1) Figurehandler.Create = Figurehandler.Create("savePath","01_Figures\"), 
    end
    if isscalar(dbstack), close all; clc; end
    f = {};
    nr_files = 500;
    
    data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\rt_cdmm_pecs\", nr_files);
    mu = groupsummary(data, ["nr_samples", "poly_order"], "mean", ["timeCDMM", "timePECS"]);
    
    poly_order_vec = sort(unique(mu.poly_order), "descend") - 1;
    nr_samples_vec = sort(unique(mu.nr_samples), "ascend");
    
    line_styles = [":", "--", "-"];
    marker_styles = ["s", "^"];
    
    f{end+1} = figure("Name","Runtime Samples");
    for i = 1:numel(poly_order_vec)
        poly_order = poly_order_vec(i);
        h(1,i) = loglog( ...
            nr_samples_vec, ...
            mu.mean_timeCDMM(mu.poly_order == poly_order+1), ...
            line_styles(i) + marker_styles(1), ...
            "DisplayName", "BSUM $Q = "+num2str(poly_order) + "$"); hold on; 
        h(2,i) = loglog( ...
            nr_samples_vec, ...
            mu.mean_timePECS(mu.poly_order == poly_order+1), ...
            line_styles(i) + marker_styles(2), ...
            "DisplayName", "PECS $Q = "+num2str(poly_order) + "$"); hold on;
    end
    % ylim([-35, 0]);
    box on; % grid on
    xticks(2.^(log2(nr_samples_vec(1)):log2(nr_samples_vec(end))))
    lgd = legend(reshape(h.', [], 1), "Location","northwest", "NumColumns", 2, "EdgeColor","white");
    lgd.set("IconColumnWidth", 15);
    xlabel("Sequence Length (-)");
    ylabel("Runtime (s)")
    ylim([1e-3, 1e3])
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'off';
    % ax.XGrid = 'on';
    % ax.XMinorGrid = 'off';
    
    
    %% Formating
    fh.PrepareSvg(f{end}, "colorList",["#D37611", "#005287"], "width", 31.56/2);
    
    if isscalar(dbstack)
        fh.GenerateOutputs();
        keyboard;
    end
end