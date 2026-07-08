close all; clearvars; clc; 

fh = Figurehandler.Create("savePath","01_Figures\");

plt_cd_improvement(fh);

plt_perf_cdmm_pecs(fh);
plt_t2q_cdmm_pecs(fh);
plt_rt_cdmm_pecs(fh);
plt_sen_cdmm(fh);
fh.GenerateOutputs();