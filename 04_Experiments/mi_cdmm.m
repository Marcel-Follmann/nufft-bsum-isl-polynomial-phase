close all; clearvars; clc; 

addpath(".\");
addpath("02_Data\");
addpath("99_Utilities\");

fh = Figurehandler.Create("savePath","01_Figures\");
f = {};
colorlist = flip(["#00934B", "#D37611", "#005287"]);

MC = 150;
nr_samples = [4;16;64]; 
nr_params = 3:6; 
max_mm_iters = 1e3; 
tol = 1e-3; 
nr_candidates = 1000;


T = table();
for k = 1:numel(nr_samples)
    for i = 1:numel(nr_params)
        params = Algorithm.Params(...
            "nr_samples", nr_samples(k), ...
            "nr_params",nr_params(i), ...
            "max_mm_iters",max_mm_iters, ...
            "tol", tol, ...
            "nr_candidates",nr_candidates ...
            );
        
        cdmm = Algorithm.CDMM.General(params);
        
        for mc = 1:MC
            [~, a0] = params.get_random_sequence();
            a(:,mc) = cdmm.run(a0);
            x(:,mc) = params.generate_signal(a(:,mc));
        end
        
        poly_order = nr_params(i)-1;
        samples = sprintf("N = %d",nr_samples(k));
        for mc = 1:MC
            for j = mc+1:MC
                y = xcorr(x(:,mc), x(:,j))/(nr_samples(k));
                % plot(20*log10(abs(y))); hold on; 
                mi = max((abs(y)));
                T = [T; table(mi, poly_order, samples)];
            end
        end
        
        
        clear a x
        fprintf("Order: %d, Length: %d, Done.\n", nr_params(i), nr_samples(k));
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
ds = (nr_params(2) - nr_params(1))/2;
xticks(nr_params-1)
xlim([nr_params(1) - ds, nr_params(end) + ds]-1);

for i = 1:numel(f)
    fh.PreparePdf(f{i});
end
ax.XMinorTick = 'off';
ax.YGrid ="on";
fh.GenerateOutputs()