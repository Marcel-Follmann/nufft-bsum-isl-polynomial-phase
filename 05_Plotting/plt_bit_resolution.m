function plt_bit_resolution(fh)
    arguments
        fh (1,1) Figurehandler.Create = Figurehandler.Create("savePath","01_Figures\"), 
    end
    if isscalar(dbstack), close all; clc; end
    f = {};
    colorlist = flip(["#00934B", "#D37611", "#005287"]);
    
    data = Datalogger.ChunkBased.logsAsFlatTable("02_Data\bit_res_3\");
    load("02_Data\bit_res_3\settings.mat");
    
    f{end+1} = figure("Name", "Bit Reolution");
    bits = 2:48;
    T = table();
    T0 = table();
    for i = 1:height(data)
        a0 = data(i,:).a .';
        x0 = params.generate_signal(a0);
        isl0 = params.compute_isl_dB(x0);
        for i_bit = 1:numel(bits)

            a = quantize(a0, bits(i_bit));
            x = params.generate_signal(a);
            T0.bit = bits(i_bit);
            T0.isl = params.compute_isl_dB(x);
            T0.isl0 = isl0;
            T = [T;T0];
        end
    end
    % group = categorical(data.nr_samples);
    err_db = pow2db(abs(db2pow(T.isl) - db2pow(T.isl0)));
    b = boxchart(T.bit, err_db + 0*T.isl - 0*T.isl0, "Notch","off", "GroupByColor",categorical(repmat("N = " + num2str(params.nr_samples) + ", Q = " + num2str(params.nr_params-1), height(T),1)));
    ylabel("Delta ISL (dB)")
    xlabel("Bit Resolution (-)");
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

function a_hat = quantize(a, b)
    delta = 2*pi / 2^b; 
    a_hat = delta * round(a / delta);
end