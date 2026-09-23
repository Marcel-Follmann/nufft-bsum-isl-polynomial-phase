function plt_t2q_cdmm_pecs(fh)
    arguments
        fh (1,1) Figurehandler.Create = Figurehandler.Create("savePath","01_Figures\"), 
    end
    if isscalar(dbstack), close all; clc; end
    f = {};
    nr_files = 8;%500;
    
    
    data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\t2q_cdmm_pecs_2\", nr_files);
    data.isl_pecs(data.isl_pecs == -Inf) = NaN;
    mu = groupsummary(data, ["poly_order", "nr_samples"], "mean", ["isl_cdmm", "isl_pecs", "psl_cdmm", "psl_pecs"]);
    
    poly_order_vec = sort(unique(mu.poly_order), "descend") - 1;
    nr_samples_vec = sort(unique(mu.nr_samples), "ascend");
    
    line_styles = [":", "--", "-"];
    marker_styles = ["s", "^"];
    
    f{end+1} = figure("Name","Time to Quality");
    for i = 1:numel(poly_order_vec)
        poly_order = poly_order_vec(i);
        h(1,i) = semilogx( ...
            nr_samples_vec, ...
            mu.mean_isl_cdmm(mu.poly_order == poly_order+1), ...
            line_styles(i) + marker_styles(1), ...
            "DisplayName", "BSUM $Q = "+num2str(poly_order) + "$"); hold on; 
        h(2,i) = semilogx( ...
            nr_samples_vec, ...
            mu.mean_isl_pecs(mu.poly_order == poly_order+1), ...
            line_styles(i) + marker_styles(2), ...
            "DisplayName", "PECS $Q = "+num2str(poly_order) + "$"); hold on;
    end
    
    
    % ylim([-35, 0]);
    box on; % grid on
    xticks(2.^(log2(nr_samples_vec(1)):log2(nr_samples_vec(end))))
    lgd = legend(reshape(h.', [], 1), "Location","northeast", "NumColumns", 2, "EdgeColor","white");
    lgd.set("IconColumnWidth", 15);
    semilogx(nr_samples_vec, -30*log10(nr_samples_vec) - 10*log10(2), "DisplayName","Limit");
    xlabel("Sequence Length (-)");
    ylabel("ISL (dB)")
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'off';
    % ax.XTick = 'off';
    % ax.XGrid = 'on';
    % ax.XMinorGrid = 'off';
    
    
    %% Formating
    fh.PrepareSvg(f{1}, "colorList",["#D37611", "#005287"], "width", 31.56/2);
    
    if isscalar(dbstack)
        fh.GenerateOutputs();
       keyboard;
    end
end