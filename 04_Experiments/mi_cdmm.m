function mi_cdmm(nr_monte_carlo, nr_samples, vec_poly_order, nr_max_iters, tolerance)
    arguments
        nr_monte_carlo (1,1) double {mustBeInteger} = 150, 
        nr_samples (:,1) double {mustBeInteger} = [4;16;64], 
        vec_poly_order (:,1) double {mustBeInteger} = 3:6,
        nr_max_iters (1,1) double {mustBeInteger} = 1e3,
        tolerance (1,1) double = 1e-3
    end
    if isscalar(dbstack), close all; clc; end
    fprintf("Mutual Interference Evaluation of CDMM\n");

    fh = Figurehandler.Create("savePath","01_Figures\");
    f = {};
    colorlist = flip(["#00934B", "#D37611", "#005287"]);
        
    T = table();
    for k = 1:numel(nr_samples)
        for i = 1:numel(vec_poly_order)
            params = Algorithm.Params(...
                "nr_samples", nr_samples(k), ...
                "nr_params",vec_poly_order(i), ...
                "max_mm_iters",nr_max_iters, ...
                "tol", tolerance ...
                );
            
            cdmm = Algorithm.CDMM.General(params);
            
            for mc = 1:nr_monte_carlo
                [~, a0] = params.get_random_sequence();
                a(:,mc) = cdmm.run(a0);
                x(:,mc) = params.generate_signal(a(:,mc));
            end
            
            poly_order = vec_poly_order(i)-1;
            samples = sprintf("N = %d",nr_samples(k));
            for mc = 1:nr_monte_carlo
                for j = mc+1:nr_monte_carlo
                    y = xcorr(x(:,mc), x(:,j))/(nr_samples(k));
                    % plot(20*log10(abs(y))); hold on; 
                    mi = max((abs(y)));
                    T = [T; table(mi, poly_order, samples)];
                end
            end
            
            
            clear a x
            fprintf("Order: %d, Length: %d, Done.\n", vec_poly_order(i), nr_samples(k));
        end
    end
    f{end+1} = figure("Name","Mutual Interference");
    group = categorical(T.samples);
    group = reordercats(group, {'N = 4','N = 16','N = 64'});
    b = boxchart(T.poly_order, 20*log10(T.mi), 'GroupByColor',group, 'Notch','on');
    ylabel("Normed Max. MI (dB)");
    xlabel("Polynomial Order (-)");
    lgd = legend("Location","best", "Box","off", "NumColumns",3);
    lgd.IconColumnWidth = 10;
    for i = 1:numel(b)
        b(i).JitterOutliers = 'on';
        b(i).MarkerStyle = '.';
        b(i).BoxFaceColor = hex2rgb(colorlist(i));
        b(i).MarkerColor = hex2rgb(colorlist(i));
    end
    ax = gca;
     ylim([-25, 0])
    ds = (vec_poly_order(2) - vec_poly_order(1))/2;
    xticks(vec_poly_order-1)
    xlim([vec_poly_order(1) - ds, vec_poly_order(end) + ds]-1);
    
    for i = 1:numel(f)
        fh.PreparePdf(f{i});
    end
    ax.XMinorTick = 'off';
    ax.YGrid ="on";
    fh.GenerateOutputs()

    if isscalar(dbstack), keyboard, end
end
