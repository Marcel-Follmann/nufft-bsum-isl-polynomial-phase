close all; clearvars; clc; 


addpath(".\");
addpath("01_Figures\");
addpath("02_Data\");
addpath("99_Utilities\utils\");

figure_handler = FigureHandler("savePath","01_Figures\");


fig = figure("Name", "NUFFT Type 2");

% Main axes
ax1 = axes;

x = linspace(-2,2,100);
n = linspace(-1.75,1.75, 10);
dn = n(2) - n(1);
f = 1;

f = @(x) -.5*cos(2*pi*f*x) + x.^2 + .5;
y = f(x);
yn = f(n);

plot(ax1, x, y); hold on; 
plot(ax1, n, yn, 'o');
axis(ax1, [-2 2 0, 5])
ax1.Box = 'on';
xticks([]);
yticks([]);
axis off

% Desired data-space box for ax2
x2 = 2*dn*[-1 1];
y2 = [2 4.5];

% Get ax1 position in figure coordinates
pos1 = ax1.Position;   % [x y w h] in normalized units
xl = ax1.XLim;
yl = ax1.YLim;

% Convert data coords → normalized figure coords
x_norm = pos1(1) + pos1(3)*(x2(1)-xl(1))/diff(xl);
y_norm = pos1(2) + pos1(4)*(y2(1)-yl(1))/diff(yl);
w_norm = pos1(3)*diff(x2)/diff(xl);
h_norm = pos1(4)*diff(y2)/diff(yl);


x_gauss = linspace(-1.75*dn,1.75*dn, 100);
n_gauss = linspace(-1.5*dn,1.5*dn, 4);
f_gauss = @(x) exp(-x.^2);
% Create inset axes
ax2 = axes('Position',[x_norm y_norm w_norm h_norm]);
plot(ax2, x_gauss, f_gauss(x_gauss)); hold on; 
plot(ax2, n_gauss, f_gauss(n_gauss), 'o');
xlim(2*dn*[-1,1])
box(ax2,'on')
xticks([]);
yticks([]);
figure_handler.PrepareEmf(fig);
% figure_handler.GenerateOutputs();