% demo_display_results_high_perfusion.m
% Written by Namgyun Lee
% Email: ggang56@gmail.com
% Started: 09/02/2018, Last modified: 10/08/2018
% Modified by NGL 06/28/2019

%% Clean slate
close all; clear all; clc;

%% Add paths
addpath(genpath('D:\ASL_project\mfiles_nam\Bloch_Flow_MT'));

%%
FontSize = 12;
color_order1 = cbrewer('qual', 'Set1', 8);
figure('Color', 'w', 'Position', [520 208 1019 590]);

%[520 208 1019 590]
%[520 173 720 625]

%% PASL
load('PASL_high_parameter_sweep');
dt_range = dt_range.';

parameter_mean = mean(nrmse*1e2, 1).';
parameter_stdev = sqrt(var(nrmse*1e2, 0, 1)).';
color = color_order1(1,:);

nrmse_range = [0 7];
maxdev_range = [0 0.04];
comptime_range = [0 3];

subplot(2,3,1);
hold on; grid on; set(gca, 'Box', 'On');
plus_stdev  = parameter_mean + parameter_stdev;
minus_stdev = parameter_mean - parameter_stdev;
hp = patch([dt_range; dt_range(end:-1:1); dt_range(1)]*1e3, [minus_stdev; plus_stdev(end:-1:1); minus_stdev(1)], color);
plot(dt_range*1e3, parameter_mean, 'Color', color, 'LineWidth', 1);
set(hp, 'edgecolor', 'none', 'FaceAlpha', 0.5);
xlabel('Timestep (msec)');
ylabel('NRMSE \pm SD (%)');
title('PASL: NRMSE');
xlim([dt_range(1) dt_range(end)]*1e3);
ylim(nrmse_range);

parameter_mean = mean(max_deviation*1e2, 1).';
parameter_stdev = sqrt(var(max_deviation*1e2, 0, 1)).';
color = color_order1(3,:);

subplot(2,3,2);
hold on; grid on; set(gca, 'Box', 'On');
plot(dt_range*1e3, parameter_mean, 'Color', color, 'LineWidth', 1);
plus_stdev  = parameter_mean + parameter_stdev;
minus_stdev = parameter_mean - parameter_stdev;
hp = patch([dt_range; dt_range(end:-1:1); dt_range(1)]*1e3, [minus_stdev; plus_stdev(end:-1:1); minus_stdev(1)], color);
set(hp, 'edgecolor', 'none', 'FaceAlpha', 0.5);
xlabel('Timestep (msec)');
ylabel('Max Deviation \pm SD (%)');
title('PASL: Max Deviation');
xlim([dt_range(1) dt_range(end)]*1e3);
ylim(maxdev_range);

parameter_mean  = mean(computation_time*1e3, 1).';
parameter_stdev = sqrt(var(computation_time*1e3, 0, 1)).';
color = color_order1(2,:);

subplot(2,3,3);
hold on; grid on; set(gca, 'Box', 'On');
plot(dt_range*1e3, parameter_mean, 'Color', color, 'LineWidth', 1);
plus_stdev  = parameter_mean + parameter_stdev;
minus_stdev = parameter_mean - parameter_stdev;
hp = patch([dt_range; dt_range(end:-1:1); dt_range(1)]*1e3, [minus_stdev; plus_stdev(end:-1:1); minus_stdev(1)], color);
set(hp, 'edgecolor', 'none', 'FaceAlpha', 0.5);
xlabel('Timestep (msec)');
ylabel('Computation Time \pm SD (msec)');
title('PASL: Computation Time');
xlim([dt_range(1) dt_range(end)]*1e3);
ylim(comptime_range);

%% PCASL
load('PCASL_high_parameter_sweep');
dt_range = dt_range.';

parameter_mean = mean(nrmse*1e2, 1).';
parameter_stdev = sqrt(var(nrmse*1e2, 0, 1)).';
color = color_order1(1,:);

subplot(2,3,4);
hold on; grid on; set(gca, 'Box', 'On');
plus_stdev  = parameter_mean + parameter_stdev;
minus_stdev = parameter_mean - parameter_stdev;
hp = patch([dt_range; dt_range(end:-1:1); dt_range(1)]*1e3, [minus_stdev; plus_stdev(end:-1:1); minus_stdev(1)], color);
plot(dt_range*1e3, parameter_mean, 'Color', color, 'LineWidth', 1);
set(hp, 'edgecolor', 'none', 'FaceAlpha', 0.5);
xlabel('Timestep (msec)');
ylabel('NRMSE \pm SD (%)');
title('PCASL: NRMSE');
xlim([dt_range(1) dt_range(end)]*1e3);
ylim(nrmse_range);

parameter_mean = mean(max_deviation*1e2, 1).';
parameter_stdev = sqrt(var(max_deviation*1e2, 0, 1)).';
color = color_order1(3,:);

subplot(2,3,5);
hold on; grid on; set(gca, 'Box', 'On');
plot(dt_range*1e3, parameter_mean, 'Color', color, 'LineWidth', 1);
plus_stdev  = parameter_mean + parameter_stdev;
minus_stdev = parameter_mean - parameter_stdev;
hp = patch([dt_range; dt_range(end:-1:1); dt_range(1)]*1e3, [minus_stdev; plus_stdev(end:-1:1); minus_stdev(1)], color);
set(hp, 'edgecolor', 'none', 'FaceAlpha', 0.5);
xlabel('Timestep (msec)');
ylabel('Max Deviation \pm SD (%)');
title('PCASL: Max Deviation');
xlim([dt_range(1) dt_range(end)]*1e3);
ylim(maxdev_range);

parameter_mean  = mean(computation_time*1e3, 1).';
parameter_stdev = sqrt(var(computation_time*1e3, 0, 1)).';
color = color_order1(2,:);

subplot(2,3,6);
hold on; grid on; set(gca, 'Box', 'On');
plot(dt_range*1e3, parameter_mean, 'Color', color, 'LineWidth', 1);
plus_stdev  = parameter_mean + parameter_stdev;
minus_stdev = parameter_mean - parameter_stdev;
hp = patch([dt_range; dt_range(end:-1:1); dt_range(1)]*1e3, [minus_stdev; plus_stdev(end:-1:1); minus_stdev(1)], color);
set(hp, 'edgecolor', 'none', 'FaceAlpha', 0.5);
xlabel('Timestep (msec)');
ylabel('Computation Time \pm SD (msec)');
title('PCASL: Computation Time');
xlim([dt_range(1) dt_range(end)]*1e3);
ylim(comptime_range);

text(-4.5 + -65.3*2, 7.6/2, '(a)', 'FontWeight', 'bold', 'FontSize', FontSize);
text(-4.5 + -65.3*2, -0.75/2, '(b)', 'FontWeight', 'bold', 'FontSize', FontSize);
text(-69.8         , 7.6/2, '(c)', 'FontWeight', 'bold', 'FontSize', FontSize);
text(-69.8         , -0.75/2, '(d)', 'FontWeight', 'bold', 'FontSize', FontSize);
text(-4.5          , 7.6/2, '(e)', 'FontWeight', 'bold', 'FontSize', FontSize);
text(-4.5          , -0.75/2, '(f)', 'FontWeight', 'bold', 'FontSize', FontSize);
%%
export_fig('../../figs/figure3', '-r864', '-tif');

