function plt_sen_cdmm(fh)
    arguments
        fh (1,1) Figurehandler.Create = Figurehandler.Create("savePath","01_Figures\"), 
    end
    if isscalar(dbstack), close all; clc; end

    f = {};
    load("02_Data\sen.mat");
    
    a = linspace(-pi/2,pi/2,1e3);
    
    for i = 1:numel(isl)
        r(i,:) = 20*log10(abs(isl(i).interpolate(a)));
    end
    
    f{end+1} = figure("Name", "Sensitivity");
    % imagesc(a/pi, samples, r)
    % shading interp
    imagesc(a/pi, samples, r);
    axis xy;
    % view(0,90)
    map = CreateColormap((["#005287", "#FFFFFF", "#D37611"]));
    colormap(map)
    c = colorbar;
    xlabel("$a_2 / \pi$ ($-$)");
    ylabel("Sequence Length ($-$)");
    ylabel(c, "ISL (dB)");
    ticks = yticks;
    ticks = [2, ticks(ticks > 2)];
    yticks(ticks);
    ylim([ticks(1), ticks(end)]);
    
    for i = 1:numel(f)
        fh.PreparePdf(f{i});
    end


    if isscalar(dbstack)
        fh.GenerateOutputs();
        keyboard;
    end
end

function map = CreateColormap(colorlist)
    N = 256;

    for iColor = 1:numel(colorlist)
        rgbList(iColor, :) = hex2rgb(colorlist(iColor));
    end

    map = interp1(linspace(0,1,iColor), rgbList, linspace(0,1,N));
end