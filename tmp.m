close all; clearvars; clc; 

addpath("99_Utilities\utils\");
addpath("01_Figures\");
figureHandler = FigureHandler("savePath","01_Figures\");
fig = {}; 

A1 = 3;
A2 = .25;
f1 = 1;
f2 = 10;
phi1 = -2*pi*0.065;
phi2 = 2*pi*0.3;
t = linspace(0,1,1e4);
x = A1*sin(2*pi*f1*t+phi1) + A1*1.5;
y = x + A2*sin(2*pi*f2*t + phi2);
xub = x+A2; 
xlb = x-A2; 

X = 2*pi*f1*A1*cos(2*pi*t+phi1);
Y = X+2*pi*A2*cos(2*pi*f2*t + phi2);
Xub = X+2*pi*A2;
Xlb = X-2*pi*A2;

tYlb = t(Xlb <= 0);
Ylb = Xlb(Xlb <= 0);
tYub = t(Xub<=0);
Yub = Xub(Xub <= 0);


fig{end+1} = figure("Name","Interpolation Effects"); hold on; 
plot(t, x, "DisplayName","Function"); 
plot(t, y, "DisplayName","Interpolation");
plot(t, xub, 'k--', "DisplayName","Error Bounds");
plot(t, xlb, 'k--', "HandleVisibility","off");
lgd = legend("Location","southwest");
lgd.ItemTokenSize = 15;
grid off; box on; 

fig{end+1} = figure("Name","Interpolation Effects 2"); hold on; 
plot(t, X, "DisplayName","Derivative"); 
plot(t, Y, "DisplayName","Interpolation");
plot(t, Xub, 'k--', "DisplayName","Error Bounds");
plot(t, Xlb, 'k--', "HandleVisibility","off");
% plot(tYub, Yub, 'r', "DisplayName","MM region");
% plot(tYlb, Ylb, 'r', "HandleVisibility","off");
% plot([t(1), t(end)], [0,0], 'k', "HandleVisibility","off")
lgd = legend("Location","southwest");
lgd.ItemTokenSize = 15;
grid on; box on; 

for i = 1:numel(fig)
    figureHandler.PreparePdf(fig{i}, "colorList",["#000000", "#000000", "#005287", "#D37611", "#00934B"]);
end
figureHandler.GenerateOutputs();
