close all; 
% clearvars; 
clc; 
% 
% data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\rt_cdmm_pecs_2\");
% settings = load("02_Data\rt_cdmm_pecs_2\settings.mat");
% fh = Figurehandler.Create();
fig = {};
% colorlist = flip(["#8E2A6B", "#00934B", "#D37611", "#005287"]);
% 
% 
% maxTime = 2;
% 
% algorithms = ["direct", "fft"];
% setup.direct = @(p) Algorithm.CDMM.Direct(p);
% setup.fft = @(p) Algorithm.CDMM.FFT(p);
% cost.direct = @(N,P) (P-1).*N.^3;                 
% cost.fft    = @(N,P) N.^(P-1);
% ref = struct("direct", [NaN NaN], "fft", [NaN NaN]);
% 
% for alg = algorithms
%     data.("timeCDMM_" + alg) = nan(height(data), 1);
%     data.("islCDMM_" + alg) = nan(height(data), 1);
%     data.("predCDMM_" + alg) = nan(height(data), 1);
% end
% 
% [groups, samples, orders] = findgroups(data.nr_samples, data.poly_order); 
% 
% for alg = algorithms
%     [~, gOrder] = sort(cost.(alg)(samples, orders));
%     for g = gOrder.'
%         idxs = find(groups == g);
% 
%         nr_samples = samples(g);
%         poly_order = orders(g);
% 
%         params = Algorithm.Params( ...
%             "nr_samples", nr_samples, ...
%             "nr_params", poly_order ...
%             );
% 
%         fprintf("Nr. samples: %d, Order: %d\n", nr_samples, poly_order-1);
% 
%         f = cost.(alg)(nr_samples, poly_order);
%         tPred = ref.(alg)(2) * f / ref.(alg)(1);
%         data.("predCDMM_" + alg)(idxs) = tPred;
% 
%         if tPred > maxTime
%             fprintf("  %-6s SKIP   (pred. %8.1f s)\n", alg, tPred);
%             continue
%         end
% 
%         instance = setup.(alg)(params);
% 
%         for i = 1:3%:numel(idxs)
%             r = idxs(i);
%             a0 = data.a0{r};
% 
%             tic; 
%             a = instance.run(a0);
%             data.("timeCDMM_" + alg)(r) = toc;
% 
%             x = params.generate_signal(a);
%             data.("islCDMM_" + alg)(r) = params.compute_isl_dB(x);
%         end
% 
% 
%         t = mean(data.("timeCDMM_" + alg)(idxs), "all", "omitmissing");
%         if isnan(ref.(alg)(1)) || f >= ref.(alg)(1)
%             ref.(alg) = [f, t];                    % extrapolate from the biggest run
%         end
%         fprintf("  %-6s %8.3f s (pred. %8.3f s)\n", alg, t, tPred);
%     end   
% end
S = groupsummary(data, ["nr_samples", "poly_order"], {"mean", "mean", "mean", "mean", "mean", "mean"}, ...
                 ["timeCDMM", "timePECS", "timeCDMM_direct", "islCDMM_direct", "timeCDMM_fft", "islCDMM_fft"]);

% --- Plot: x = nr_samples, eine Linie pro poly_order ---


for p = unique(S.poly_order).'
    T = sortrows(S(S.poly_order == p, :), "nr_samples");
    fig{end+1} = figure("Name", "Runtime all"); 
    for alg = algorithms
        name = sprintf("%s", alg);
        loglog(T.nr_samples, T.("mean_timeCDMM_" + alg), "-o", "DisplayName", name); hold on;
    end
    name = sprintf("general");
    loglog(T.nr_samples, T.("mean_timeCDMM"), "-o", "DisplayName", name); hold on;
    name = sprintf("PECS");
    loglog(T.nr_samples, T.("mean_timePECS"), "-o", "DisplayName", name); hold on;

    title(sprintf("Order %d", p-1));
    nr_samples_vec = unique(T.nr_samples);
    xticks(2.^(log2(nr_samples_vec(1)):log2(nr_samples_vec(end))))
    xlabel("Nr. samples"); ylabel("Mean runtime [s]");
    legend("Location", "best"); 
    grid("on"); 
end



for p = unique(S.poly_order).'
    T = sortrows(S(S.poly_order == p, :), "nr_samples");
    fig{end+1} = figure("Name", "ISL direct fft");
    for alg = algorithms
        name = sprintf("%s", alg);
        loglog(T.nr_samples, T.("mean_islCDMM_" + alg), "-o", "DisplayName", name); hold on;
    end

    title(sprintf("Order %d", p-1));
    nr_samples_vec = unique(T.nr_samples);
    xticks(2.^(log2(nr_samples_vec(1)):log2(nr_samples_vec(end))))
    xlabel("Nr. samples"); ylabel("Mean runtime [s]");
    legend("Location", "best"); 
    grid("on"); 
end


for samps = [16, 32, 64]
    fig{end+1} = figure("Name", "Order wise comparison");
    T = sortrows(S(S.nr_samples == samps, :), "poly_order");
    for alg = algorithms
        name = sprintf("%s", alg);
        semilogy(T.poly_order-1, T.("mean_timeCDMM_" + alg), "-o", "DisplayName", name); hold on;
    end
    name = sprintf("general");
    semilogy(T.poly_order-1, T.("mean_timeCDMM"), "-o", "DisplayName", name); hold on;
    name = sprintf("PECS");
    semilogy(T.poly_order-1, T.("mean_timePECS"), "-o", "DisplayName", name); hold on;
    
    title(sprintf("Nr. samples %d", samps));
    poly_order_vec = unique(T.poly_order);
    xticks(poly_order_vec(1):poly_order_vec(end))
    xlabel("Polynomial Order"); ylabel("Mean runtime [s]");
    legend("Location", "best"); 
    grid("on"); 
end


for i = 1:3
    fh.PrepareSvg(fig{i}, "colorList",colorlist)
    fh.PrepareSvg(fig{i+3}, "colorList",colorlist(1:2))
    fh.PrepareSvg(fig{i+6}, "colorList",colorlist)
end




function completeColorList = CreateColorList(incompleteColorList, listLength)
	% Calculate number of complete iterates
	paletteLength = numel(incompleteColorList);
	shadingIterates = ceil(listLength / paletteLength);

	% Calculate shading step size
	maximalShadingFactor = .75; % 0: no shading | 1: full shade
	shadingStepFactor = maximalShadingFactor / shadingIterates;
	
	% Iterate over palette
	completeColorList = strings(1,listLength);
	for i = 1:listLength
		iteration = ceil(i/paletteLength);
		paletteIndex = mod(i-1,paletteLength)+1;
		baseColor = Figurehandler.hex2rgb(incompleteColorList(paletteIndex));

		% Increase shading with each new iterate
		shading = shadingStepFactor * (iteration-1);
		rgbShaded = baseColor * (1 - shading);
		completeColorList(i) = Figurehandler.rgb2hex(rgbShaded);
	end
end
