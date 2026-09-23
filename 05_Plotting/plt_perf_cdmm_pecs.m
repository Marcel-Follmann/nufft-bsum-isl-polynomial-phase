function plt_perf_cdmm_pecs(fh)
    arguments
        fh (1,1) Figurehandler.Create = Figurehandler.Create("savePath","01_Figures\"), 
    end
    if isscalar(dbstack), close all; clc; end
    fig = {};
    
    cdmm = load("02_Data\perf_2\perf_cdmm.mat");
    pecs = load("02_Data\perf_2\perf_pecs.mat");
    
    % Reconstruct a over the Optimization
    [a, a_mm] = reorder_cdmm(cdmm);
    params = pecs.params;
    
    x_cdmm = params.generate_signal(a);
    x_pecs = params.generate_signal(pecs.a);
    isl_dB_cdmm = params.compute_isl_dB(x_cdmm);
    isl_dB_pecs = params.compute_isl_dB(x_pecs);
    
    % a_mm = mod(a_mm, pi);
    
    fig{end+1} = figure("Name", "ISL in Coordinate"); 
    % theta = linspace(0, 2*pi, 1e4);
    p = numel(cdmm.isl)+1;
    theta = linspace(cdmm.isl(end).s_range(1), cdmm.isl(end).s_range(2), 1e6);
    plot(a_mm(end,:)/pi, isl_dB_cdmm(end,:), 'x', "DisplayName","Minima"); hold on; 
    plot(theta/(pi), 10*log10(real(cdmm.isl(end).interpolate(theta))), '-', "DisplayName","ISL");hold on; 
    xlim(a_mm(end,:)/pi + 1e-2*[-1, 1])
    xlabel("$a_" + num2str(p) + " / \pi$ (-)");
    ylabel("ISL (dB)");
    lgd = legend("Location","northwest", "Box", "off", "NumColumns", 2);
    lgd.IconColumnWidth = 10;
    box on; %grid on; 
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'off';
    
    fig{end+1} = figure("Name", "Convergence Plot"); 
    semilogx(isl_dB_pecs, '-', "DisplayName","PECS"); hold on;
    semilogx(isl_dB_cdmm, '-', "DisplayName","BSUM");
    box on; %grid on; 
    ylabel("ISL (dB)");
    xlabel("Number of MM steps (-)")
    lgd = legend("Location","northeast", "EdgeColor","white", "NumColumns", 1);
    lgd.IconColumnWidth = 10;
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'off';
    
    % fig{end+1} = figure("Name","ACF");
    % N = cdmm.params.nr_samples;
    % f = -N+1:N-1;
    % 
    % 
    
    % r_start = params.compute_autocorr(x_pecs(:,1) ,"symmetric", true);
    % r_cdmm = params.compute_autocorr(x_cdmm(:,end) ,"symmetric", true);
    % r_pecs = params.compute_autocorr(x_pecs(:,end) ,"symmetric", true);
    % 
    % plot(f,20*log10(abs(r_start)), "DisplayName","Initial"); hold on; 
    % plot(f,20*log10(abs(r_cdmm)), "DisplayName","BSUM"); 
    % plot(f,20*log10(abs(r_pecs)), "DisplayName","PECS");
    % ylim([-65, 5]);
    % ylabel("Magnitude (dB)");
    % xlabel("Lags (-)");
    % lgd = legend("Location","northeast", "EdgeColor","white", "NumColumns", 1);
    % lgd.IconColumnWidth = 10;
    % ax = gca;
    % ax.YGrid = 'on';
    % ax.YMinorGrid = 'off';
    
    fh.PrepareSvg(fig{1}, "markerSize", 5, "width", 31.56/2);%, "colorList",["#005287", "#005287", "#005287", "#005287", "#D37611", "#005287", "#00934B", "#CCCCCC"])
    for i = 2:numel(fig)
        fh.PrepareSvg(fig{i}, "markerSize", 5, "width", 31.56/2);%, "colorList",flip(["#005287", "#D37611", "#00D16C", "#00934B", "#00753D", "#004725"]));
    end
    clear ans
    if isscalar(dbstack)
        fh.GenerateOutputs();
        keyboard;
    end
end

function [a, a_mm] = reorder_cdmm(cdmm)
    a = cdmm.a(:,end-1);
    a_mm = cdmm.a_MM(:, find(cdmm.i_MM == 0,1,"last"):end);
    % a_mm = a_mm(:,:,end,end);
    a_mm = permute(a_mm, [3, 2, 1]);
    a = a.*ones(cdmm.params.nr_params, size(a_mm,2), 1, cdmm.params.nr_cd_iterations);
    a(end,:,:) = a_mm;
    
    a_mm = permute(a_mm, [2, 3, 1]);
end