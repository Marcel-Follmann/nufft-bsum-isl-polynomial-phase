function plt_cd_improvement(fh)
    arguments
        fh (1,1) Figurehandler.Create = Figurehandler.Create("savePath","01_Figures\"), 
    end
    if isscalar(dbstack), close all; clc; end
    f = {};
    colorlist = flip(["#00934B", "#D37611", "#005287"]);
    
    data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\ci\");
    load("02_Data\ci\settings.mat");
    
    f{end+1} = figure("Name", "ISL Improvement");
    group = categorical(data.nr_samples);
    % group = reordercats(group, {'N = 16','N = 32', 'N = 64'});
    b = boxchart(data.coordinate, -data.improvement, "GroupByColor",data.nr_samples, "Notch","off");
    ylabel("ISL Improvement (dB)")
    xlabel("Coordinate (-)");
    
    ax = gca;
    ds = 1/2;
    xticks(2:poly_order-1)
    xlim([2 - ds, poly_order + ds-1]);
    
    for i = 1:numel(b)
        b(i).JitterOutliers = 'on';
        b(i).MarkerStyle = '.';
        b(i).BoxFaceColor = hex2rgb(colorlist(i));
        b(i).MarkerColor = hex2rgb(colorlist(i));
    end
    box on;
    lgd = legend("Location","best", "EdgeColor","white", "NumColumns", 1);
    lgd.set("IconColumnWidth", 15);
    fh.PrepareSvg(f{end});
    ax.XMinorTick = 'off';
    ax.YGrid ="on";

    if isscalar(dbstack)
        fh.GenerateOutputs();
        keyboard;
    end
end