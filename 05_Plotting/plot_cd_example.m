close all; clearvars; clc; 

addpath(".\");
addpath("01_Figures\");
addpath("02_Data\");
addpath("99_Utilities\utils\");

figure_handler = FigureHandler("savePath","01_Figures\");
fig = {};

%% Params
f = 2;
rho = linspace(0,1,f*10);
phi = linspace(0,2*pi, f*20).';

z_polar =@(rho, phi) rho.^2.5 .* ones(size(phi)); 

x = rho .* cos(phi);
y = rho .* sin(phi);

% [X, Y, Z] = meshgrid(x, y, z);
%% Simulation
fig{end+1} = figure("Name", "CD Example");
surf(x + .5, y+ .5, z_polar(rho, phi), 'FaceColor',.95*[1,1,1], ...
    'EdgeColor','k', ...
    'FaceAlpha',.95, ...
    'LineWidth', 2);
grid off;
colormap([0,0,0])
% lighting none
% shading flat
axis off
% view(3)

for i = 1:numel(fig)
    figure_handler.PrepareEmf(fig{i});
end
figure_handler.GenerateOutputs();
%% Visualization