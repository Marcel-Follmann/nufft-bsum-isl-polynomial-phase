close all; clearvars; clc; 

addpath(".\");
addpath("02_Data\");
addpath("99_Utilities");
rng("shuffle");

fh = Figurehandler.Create("savePath","01_Figures\");
f = {};
colorlist = flip(["#00934B", "#D37611", "#005287"]);

tol = 1e-7;
MC = 100;
nr_samples = [16; 32; 64];
nr_params = 5;



T = table();
for k = 1:numel(nr_samples)
    for mc = 1:MC
        params = Algorithm.Params(...
        "nr_samples", nr_samples(k), ...
        "nr_params",nr_params, ...
        "max_mm_iters",1e4, ...
        "tol", tol, ...
        "nr_candidates",1001, ...
        "nr_cd_iterations", 1 ...
        );
        
        cdmm = Algorithm.CDMM.General(params);
        dl = Datalogger.Create("file_name","Coordinate_Improvement");
        [~, a0] = params.get_random_sequence();
        cdmm.run(a0, "data_logger", dl);
        pause(0.1);
        load("02_Data\Coordinate_Improvement.mat");
        pause(0.1);
        for i = 1:size(a,2)
            x = params.generate_signal(a(:,i));
            isl_coordinate(i) = params.compute_isl_dB(x);
        end
        tmp = diff(isl_coordinate);
        for i = 1:size(a,2)-1
            coordinate = i + 1;
            samples = sprintf("N = %d",nr_samples(k));
            improvement = tmp(i);
            T = [T; table(improvement, coordinate, samples)];
        end
    end
end

f{end+1} = figure("Name", "ISL Improvement");
group = categorical(T.samples);
group = reordercats(group, {'N = 16','N = 32', 'N = 64'});
b = boxchart(T.coordinate, -T.improvement, "GroupByColor",T.samples, "Notch","off");
ylabel("ISL Improvement (dB)")
xlabel("Coordinate (-)");

ax = gca;
ds = 1/2;
xticks(2:nr_params-1)
xlim([2 - ds, nr_params + ds-1]);

for i = 1:numel(b)
    b(i).JitterOutliers = 'on';
    b(i).MarkerStyle = '.';
    b(i).BoxFaceColor = hex2rgb(colorlist(i));
    b(i).MarkerColor = hex2rgb(colorlist(i));
end
lgd = legend("Location","best", "Box", "off");

fh.PreparePdf(f{1});
ax.XMinorTick = 'off';
ax.YGrid ="on";
fh.GenerateOutputs();


