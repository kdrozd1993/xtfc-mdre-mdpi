clc
clear all
close all

load('Results_MDRE/mdre_xtfc_mc_w-11_b01_T01_20220317.mat')
% load('Results_MDRE/mdre_xtfc_mc_w-11_b01_T001_20220314.mat');
results1 = results;
%mcn = options.amountMcSimulations;

load('Results_MDLE/mdle_xtfc_mc_w-11_b01_T01_20220317.mat')
% load('Results_MDLE/mdle_xtfc_mc_w-11_b01_T001_20220315.mat');
results2 = results;

load('Results_MDRE/mdre_xtfc_w-11_b01_T01_halton_20220401.mat')
% load('Results_MDRE/mdre_xtfc_w-11_b01_T001_halton_20220401.mat');
results3 = results;

load('Results_MDLE/mdle_xtfc_w-11_b01_T01_halton_20220401.mat');
% load('Results_MDLE/mdle_xtfc_w-11_b01_T001_halton_20220401.mat');
results4 = results;

time = results1(1).time.interp;

set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
%set(groot,'defaultLineLineWidth',2)
set(0,'DefaultAxesFontSize', 10)
%set(0,'DefaultLineMarkerSize',10);

f1 = figure(1);

int = floor(linspace(1, length(time), 30));

p1 = plot(time, results1(1).residualNorm.LyapKE, '-', 'Color', [0 100/255 0], 'LineWidth', 3); hold on
p2 = plot(time, results1(1).residualNorm.RK4KE, 'k-', 'LineWidth', 3); hold on

for i = 1:mcn
    res1(:,i) = results1(i).residualNorm.PinnKE;
    res2(:,i) = results2(i).residualNorm.PinnLKE;
    %semilogy(time, res1(:,i), 'b-'); hold on
end

[b1, tf1] = rmoutliers(res1(int(end-1),:), 'percentiles', [0, 90]);
res1(:, tf1) = [];
meanRes1 = mean(res1, 2);
stdRes1 = std(res1, 0, 2);
minRes1 = min(res1, [], 2);
maxRes1 = max(res1, [], 2);
ci951 = 1.96*stdRes1./sqrt(mcn);

[b2, tf2] = rmoutliers(res2(int(end-1),:), 'percentiles', [0, 90]);
res2(:, tf2) = [];
meanRes2 = mean(res2, 2);
stdRes2 = std(res2, 0, 2);
minRes2 = min(res2, [], 2);
maxRes2 = max(res2, [], 2);
ci952 = 1.96*stdRes2./sqrt(mcn);

p3 = plot(time, meanRes1, 'b-', 'LineWidth', 4); hold on
% fill([time(1:end-1)' fliplr(time(1:end-1)')], [meanRes1(1:end-1)'+ci951(1:end-1)', fliplr(meanRes1(1:end-1)'-ci951(1:end-1)')], 'b', 'FaceAlpha', 0.5, 'LineStyle', 'none'); hold on
% fill([time(1:end-1)' fliplr(time(1:end-1)')], [meanRes2(1:end-1)'+ci952(1:end-1)', fliplr(meanRes2(1:end-1)'-ci952(1:end-1)')], [0 100/255 0], 'FaceAlpha', 0.5, 'LineStyle', 'none'); hold on
fill([time(1:end-1)' fliplr(time(1:end-1)')], [maxRes1(1:end-1)', fliplr(minRes1(1:end-1)')], 'b', 'FaceAlpha', 0.5, 'LineStyle', 'none'); hold on


res3 = results3.residualNorm.PinnKE;
res4 = results4.residualNorm.PinnLKE;

p5 = plot(time, res3, 'c-', 'LineWidth', 2); hold on
p6 = plot(time, res4, 'y-', 'LineWidth', 2); hold on

p4 = plot(time, meanRes2, 'r-', 'LineWidth', 4); hold on
fill([time(1:end-1)' fliplr(time(1:end-1)')], [maxRes2(1:end-1)', fliplr(minRes2(1:end-1)')], 'r', 'FaceAlpha', 0.5, 'LineStyle', 'none'); hold on

%p6 = plot(time, res4, 'y-', 'LineWidth', 2); hold on
%fill([time(1:end-1)' fliplr(time(1:end-1)')], [maxRes2(1:end-1)', fliplr(minRes2(1:end-1)')], 'r', 'FaceAlpha', 0.5, 'LineStyle', 'none'); hold on

h1 = plot(NaN, NaN, '-', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2); hold off

% legend([p3, p4, p5, h1, p1, p2], 'PI-MDRE unif.', ...
%     'PI-MDLE unif.', 'PI-MDRE halt.', 'PI-MDLE halt.', 'NG', ...
%     'RK4', 'Location', 'South', 'NumColumns', 3);

legend([p3, p4, p5, h1, p1, p2], 'PI-MDRE unif.', ...
    'PI-MDLE unif.', 'PI-MDRE halt.', 'PI-MDLE halt.', 'NG', ...
    'RK4', 'Location', 'North', 'NumColumns', 2);

ylabel('Normwise Error')
xlabel('Time (s)')
% title('h = 0.001 s')
title('h = 0.01 sec')
grid
set(gca, 'YScale','log')
% set(f1, 'Position',  [100, 100, 800, 500])
set(f1, 'Units', 'inches', 'Position',  [1, 1, 5, 3])
set(gcf, 'PaperPositionMode', 'auto');