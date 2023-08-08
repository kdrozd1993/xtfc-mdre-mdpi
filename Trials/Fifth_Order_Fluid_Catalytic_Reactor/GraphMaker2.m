clc
clear all
close all

load('Results_MDRE/mdre_xtfc_vary_horizon_20220417.mat')
results1 = results;
load('Results_MDLE/mdle_xtfc_vary_horizon_20220417.mat')
results2 = results;

time = results1(1).time.interp;

set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
%set(groot,'defaultLineLineWidth',2)
set(0,'DefaultAxesFontSize', 10)
%set(0,'DefaultLineMarkerSize',10);

%mcn = options.amountMcSimulations;

f1 = figure(1);

int = floor(linspace(1, length(time), 30));

p1 = plot(time, results1(1).residualNorm.LyapKE, 'r-', 'LineWidth', 2); hold on
p2 = plot(time, results1(1).residualNorm.RK4KE, 'k-', 'LineWidth', 2); hold on

res1 = results1.residualNorm.PinnKE;
res2 = results2.residualNorm.PinnLKE;

p3 = plot(time, res1, 'c-', 'LineWidth', 2); hold on
p4 = plot(time, res2, '-', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2); hold on

legend([p3, p4, p1, p2], 'PI-MDRE halt.', 'PI-MDLE halt.', 'NG', 'RK4', 'Location', 'NorthWest');
ylabel('Normwise Error')
xlabel('Time (s)')
title('h varies from 0.001 s to 0.01 s')
grid
set(gca, 'YScale','log')
% set(f1, 'Position',  [100, 100, 800, 500]);
set(f1, 'Units', 'inches', 'Position',  [1, 1, 5, 3])
set(gcf, 'PaperPositionMode', 'auto');