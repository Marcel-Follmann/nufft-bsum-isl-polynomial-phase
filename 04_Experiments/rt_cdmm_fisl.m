function rt_cdmm_fisl(nr_monte_carlo, nr_samples)
    arguments
        nr_monte_carlo (1,1) double {mustBeInteger} = 10, 
        nr_samples (:,1) double {mustBeInteger} = 32:32:256, 
    end
    if isscalar(dbstack), close all; clc; end
    fprintf("Runtime anylysis: CDMM vs. FISL\n");

    fh = Figurehandler.Create("savePath","01_Figures\");
    f = {};
    colorlist = flip(["#00934B", "#D37611", "#005287"]);

    T = table();
    %% Simulation
    for i = 1:numel(nr_samples)
        params = Algorithm.Params(...
        "nr_samples", nr_samples(i), ...
        "nr_params", 3 ...
        );
    
        % Initialization
        cdmm = Algorithm.CDMM.General(params);
        cdmm_ood = Algorithm.CDMM.OOD(params);
        fisl = Algorithm.FISL(params);
        
        fprintf("P: %d ...", nr_samples(i));
        
        x0 = zeros(nr_samples(i),nr_monte_carlo);
        a0 = zeros(3,nr_monte_carlo);
        for mc = 1:nr_monte_carlo
            [x0(:,mc), a0(:,mc)] = params.get_random_sequence();
        end
    
        time = zeros(nr_monte_carlo, 1);
        
        fprintf(" CDMM-OOD ...");
        type = repmat("CDMM-OOD", numel(time), 1);
        samples = nr_samples(i) * ones(numel(time), 1);
        for mc = 1:nr_monte_carlo 
            tic; cdmm_ood.run(a0(:,mc), "a_range",a0(3,mc)+[-.5, .5]); time(mc) = toc; 
            time_median(mc) = median(time(1:mc));
            nexttile(1); loglog(time_median(1:mc)); title(type(1)); nexttile(2); histogram(time(1:mc), ceil(10*log10(mc)+eps)); drawnow;
        end
        T = [T;table(type, samples, time)];
    
    
        fprintf(" CDMM ...");
        type = repmat("CDMM", numel(time), 1);
        samples = nr_samples(i) * ones(numel(time), 1);
        for mc = 1:nr_monte_carlo 
            tic; cdmm.run(a0(:,mc), "a_range",a0(3,mc)+[-.5, .5]); time(mc) = toc; 
            time_median(mc) = median(time(1:mc));
            nexttile(1); loglog(time_median(1:mc)); title(type(1)); nexttile(2); histogram(time(1:mc), ceil(10*log10(mc)+eps)); drawnow;
        end
        T = [T;table(type, samples, time)];
    
    
        fprintf(" FISL ...");
        type = repmat("FISL", numel(time), 1);
        samples = nr_samples(i) * ones(numel(time), 1);
        for mc = 1:nr_monte_carlo
            tic; fisl.run(x0(:,mc)); time(mc) = toc; 
            time_median(mc) = median(time(1:mc));
            nexttile(1); loglog(time_median(1:mc)); title(type(1)); nexttile(2); histogram(time(1:mc), ceil(10*log10(mc)+eps)); drawnow;
        end
        T = [T;table(type, samples, time)];
    
    
        fprintf(" Done.\n");
    end
    close all; 
    
    %% Visualization
    f{end+1} = figure("Name","Runtime FISL vs CDMM A");
    b = boxchart(T.samples, T.time, 'GroupByColor',T.type, 'Notch','on');
    xlabel("Sequence Length (-)");
    ylabel("Runtime (s)");
    ds = (nr_samples(2) - nr_samples(1))/2;
    xticks(nr_samples)
    xlim([nr_samples(1) - ds, nr_samples(end) + ds]);
    ax = gca;
    ax.XMinorTick = 'off';
    lgd = legend("Location","northwest", "Box","off");
    lgd.IconColumnWidth = 10;
    
    for i = 1:numel(b)
        b(i).JitterOutliers = 'on';
        b(i).MarkerStyle = '.';
        b(i).BoxFaceColor = hex2rgb(colorlist(i));
        b(i).MarkerColor = hex2rgb(colorlist(i));
    end
    
    fh.PreparePdf(f{1});
    ax.YScale ="log";
    ax.XMinorTick = 'off';
    fh.GenerateOutputs();

   if isscalar(dbstack), keyboard, end
end



