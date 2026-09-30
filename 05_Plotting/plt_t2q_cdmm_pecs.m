function plt_t2q_cdmm_pecs(fh)
    arguments
        fh (1,1) Figurehandler.Create = Figurehandler.Create("savePath","01_Figures\"), 
    end
    if isscalar(dbstack), close all; clc; end
    f = {};
    
    
    data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\t2q_cdmm_pecs_2\");
    data.isl_pecs(data.isl_pecs == -Inf) = NaN;
    mu = groupsummary(data, ["poly_order", "nr_samples"], "mean", ["isl_cdmm", "isl_pecs", "isl_fisl", "psl_cdmm", "psl_pecs", "psl_fisl"]);
    
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
    semilogx(nr_samples_vec, mu.mean_isl_fisl(mu.poly_order == poly_order +1), "DisplayName","FISL");
    xlabel("Sequence Length (-)");
    ylabel("ISL (dB)")
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'off';
    % ax.XTick = 'off';
    % ax.XGrid = 'on';
    % ax.XMinorGrid = 'off';

    f{end+1} = figure("Name","PSL results");
    for i = 1:numel(poly_order_vec)
        poly_order = poly_order_vec(i);
        h(1,i) = semilogx( ...
            nr_samples_vec, ...
            mu.mean_psl_cdmm(mu.poly_order == poly_order+1), ...
            line_styles(i) + marker_styles(1), ...
            "DisplayName", "BSUM $Q = "+num2str(poly_order) + "$"); hold on; 
        h(2,i) = semilogx( ...
            nr_samples_vec, ...
            pow2db(sqrt(nr_samples_vec/2).*sqrt(db2pow(mu.mean_isl_cdmm(mu.poly_order == poly_order+1)))), ...
            line_styles(i) + marker_styles(2), ...
            "DisplayName", "Bound $Q = "+num2str(poly_order) + "$"); hold on;
    end
    % semilogx(nr_samples_vec, mu.mean_psl_fisl(mu.poly_order == poly_order +1), "DisplayName","FISL");
    box on; % grid on
    xticks(2.^(log2(nr_samples_vec(1)):log2(nr_samples_vec(end))))
    ylim([-25, 10])
    lgd = legend(reshape(h.', [], 1),"Location","northeast", "NumColumns", 2, "EdgeColor","white");
    lgd.set("IconColumnWidth", 15);
    
    xlabel("Sequence Length (-)");
    ylabel("PSL (dB)")
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'off';
    
    
    %% Formating
    fh.PreparePdf(f{1}, "colorList",["#000000", CreateColorList(["#D37611", "#005287"], 6)]);
    fh.PreparePdf(f{2}, "colorList", CreateColorList(["#000000", "#005287"], 6));
    
    if isscalar(dbstack)
        fh.GenerateOutputs();
        keyboard;
    end
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