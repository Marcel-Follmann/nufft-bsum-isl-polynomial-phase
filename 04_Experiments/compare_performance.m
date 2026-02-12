close all; clearvars; clc; 

addpath(".\");
addpath("02_Data\");
rng("shuffle");

nSamples = 20;
nPolyOrder = 3;

pecs = PECS("nr_samples", nSamples, "nr_params",nPolyOrder, "save_data", true, "file_name", "pecs_performance");
cdmm = CDMMT3("nr_samples", nSamples, "nr_params",nPolyOrder, "save_data", true, "file_name", "cdmm_performance", ...
    "max_mm_iters",1e3, "nr_candidates",5);

[x0, a0] = pecs.getRandomSequence();
a0 = mod(a0, pi/2); % as this is periodic
pecs.mpsl(x0, a0);
cdmm.type2(a0);

delete(cdmm);
delete(pecs);


