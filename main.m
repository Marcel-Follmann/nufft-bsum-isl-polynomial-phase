close all; clearvars; clc; 
addpath("99_Utilities\utils\");

figureHandler = FigureHandler("savePath","01_Figures\");
figures = {};

%% Parameters
maxIterations = 1;
N = 100;
Q = 4;

a0(:,1) = randn(Q,1);
% a0 = [randn(3,1);zeros(Q-3,1)  ]; % Tends to produce better ISLs

%% Coordinate Descent Simulation

cd_nufft_mm = CdNufftMm( ...
	"nr_samples",N,...
    "nr_params",Q ...
	);

isl_dB(1) = cd_nufft_mm.ComputeIsl_dB(a0);

tic;
a = cd_nufft_mm.simulate(a0);
toc;
isl_db(2) = cd_nufft_mm.ComputeIsl_dB(a);


%% Plot


%% Format Figures
for iFigure = 1:numel(figures)
	figureHandler.PreparePdf(figures{iFigure});
end
% figureHandler.GenerateOutputs();
clear ans
