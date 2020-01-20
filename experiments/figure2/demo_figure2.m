% demo_figure2.m
% Written by Namgyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 06/28/2019, Last modified: 11/04/2019

%% Clean slate
close all; clear all; clc;

%% Add paths
addpath(genpath('D:\ASL_project\mfiles_nam\Bloch_Flow_MT'));

%% Define simulation parameters
% Used the recommended parameters from the consensus paper
F            = 0.6;             % brain blood flow [mL/g/min] WM: 34 [ml/100 ml/min], GM: 42 [ml/100 ml/min]
T1b          = 1650e-3;         % longitudinal relaxation time of arterial blood at 3T [sec]
T1           = 1200e-3;         % spin-lattice relaxation time of human brain at 3T [sec] WM: 1110 msec, GM: 1470 msec
T2           = 1e10;            % spin-spin relaxation time of human brain at 3T [sec]
lambda       = 0.9;             % brain/blood partition coefficient [mL of blood /g of brain tissue]
TD           = 700e-3;          % transit delay or arterial transit time (ATT) (deltat) [sec]
alpha_PASL   = 0.98;            % labeling efficiency for PASL
alpha0_PASL  = 2 * alpha_PASL;  % labelling efficiency; 1 for saturation, 2 for inversion, and 0 for control
alpha_PCASL  = 0.85;            % labeling efficiency for PCASL
alpha0_PCASL = 2 * alpha_PCASL; % labelling efficiency; 1 for saturation, 2 for inversion, and 0 for control
TW_PASL      = 800e-3;          % PASL TI1, bolus duration of the labeling pulse (tau) [sec]
TW_PCASL     = 1800e-3;         % bolus duration of the labeling pulse (tau) [sec]
m0           = 1;               % equilibrium z magnetization
m0b          = m0 / lambda;     % M0 arterial blood = M0 tissue / lambda
df           = 0;               % tissue off-resonance frequency [Hz]

T            = 4;               % duration of a signal evolution [sec]
dt1          = 3e-3;            % 1st time interval [sec]
dt2          = 35e-3;           % 4th time interval [sec]

t_labeling   = 0;               % start time of labeling RF pulses [sec]

%% Calculate miscellaneous things
% Calculate the number of time intervals
N1 = floor(T / dt1);
N2 = floor(T / dt2);

% Calculate the sampling time
t1  = (0:N1-1).' * dt1; % [sec]
t2  = (0:N2-1).' * dt2; % [sec]

nr_labeling = length(t_labeling);
labeling_on = true(nr_labeling,1);

%% Calculate the pulsed ASL signal s(t)
PASL_analytic1 = calculate_GKM_pulsed_ASL_signal(t1, F, lambda, m0b, alpha_PASL, T1, T1b, t_labeling, TD, TW_PASL);
PASL_analytic2 = calculate_GKM_pulsed_ASL_signal(t2, F, lambda, m0b, alpha_PASL, T1, T1b, t_labeling, TD, TW_PASL);

%% Calculate the continuous ASL signal s(t)
PCASL_analytic1 = calculate_GKM_continuous_ASL_signal(t1, F, lambda, m0b, alpha_PCASL, T1, T1b, t_labeling, TD, TW_PCASL);
PCASL_analytic2 = calculate_GKM_continuous_ASL_signal(t2, F, lambda, m0b, alpha_PCASL, T1, T1b, t_labeling, TD, TW_PCASL);

%% Calculate the pulsed ASL bolus signal
pulsed_ASL_bolus1 = calculate_pulsed_ASL_bolus(t1, F, lambda, m0, alpha0_PASL, T1b, t_labeling, TD, TW_PASL, labeling_on); % [magnetization/g tissue/sec]
pulsed_ASL_bolus2 = calculate_pulsed_ASL_bolus(t2, F, lambda, m0, alpha0_PASL, T1b, t_labeling, TD, TW_PASL, labeling_on); % [magnetization/g tissue/sec]

%% Calculate the continuous ASL bolus signal
continuous_ASL_bolus1 = calculate_continuous_ASL_bolus(t1, F, lambda, m0, alpha0_PCASL, T1b, t_labeling, TD, TW_PCASL, labeling_on); % [magnetization/g tissue/sec]
continuous_ASL_bolus2 = calculate_continuous_ASL_bolus(t2, F, lambda, m0, alpha0_PCASL, T1b, t_labeling, TD, TW_PCASL, labeling_on); % [magnetization/g tissue/sec]

%% Perform Bloch-Flow simulation
tau1 = dt1 * ones(N1,1, 'double'); % [sec]
tau2 = dt2 * ones(N2,1, 'double'); % [sec]

b1  = zeros(N1,1, 'double'); % [G]
gr  = zeros(N1,1, 'double'); % [G/cm]
tic; [mx1,my1,mz1] = bloch_flow(b1, gr, tau1, tau1, T1, T2, 0, 0, 2, 0, 0, 1, F, lambda, m0, zeros(N1,1)); toc;
tic; [mx2,my2,mz2] = bloch_flow(b1, gr, tau1, tau1, T1, T2, 0, 0, 2, 0, 0, 1, F, lambda, m0, pulsed_ASL_bolus1); toc;
PASL_numeric1 = mz1 - mz2; % control - labeled

b1  = zeros(N2,1, 'double'); % [G]
gr  = zeros(N2,1, 'double'); % [G/cm]
tic; [mx1,my1,mz1] = bloch_flow(b1, gr, tau2, tau2, T1, T2, 0, 0, 2, 0, 0, 1, F, lambda, m0, zeros(N2,1)); toc;
tic; [mx2,my2,mz2] = bloch_flow(b1, gr, tau2, tau2, T1, T2, 0, 0, 2, 0, 0, 1, F, lambda, m0, pulsed_ASL_bolus2); toc;
PASL_numeric2 = mz1 - mz2; % control - labeled

b1  = zeros(N1,1, 'double'); % [G]
gr  = zeros(N1,1, 'double'); % [G/cm]
tic; [mx1,my1,mz1] = bloch_flow(b1, gr, tau1, tau1, T1, T2, 0, 0, 2, 0, 0, 1, F, lambda, m0, zeros(N1,1)); toc;
tic; [mx2,my2,mz2] = bloch_flow(b1, gr, tau1, tau1, T1, T2, 0, 0, 2, 0, 0, 1, F, lambda, m0, continuous_ASL_bolus1); toc;
PCASL_numeric1 = mz1 - mz2; % control - labeled

b1  = zeros(N2,1, 'double'); % [G]
gr  = zeros(N2,1, 'double'); % [G/cm]
tic; [mx1,my1,mz1] = bloch_flow(b1, gr, tau2, tau2, T1, T2, 0, 0, 2, 0, 0, 1, F, lambda, m0, zeros(N2,1)); toc;
tic; [mx2,my2,mz2] = bloch_flow(b1, gr, tau2, tau2, T1, T2, 0, 0, 2, 0, 0, 1, F, lambda, m0, continuous_ASL_bolus2); toc;
PCASL_numeric2 = mz1 - mz2; % control - labeled

%% Display results
color_order1 = [     0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840];

background_color = [255 255 255] / 255;

LineWidth = 1;
FontSize = 10;

spacing = 0.005;
vspacing = 0.01;

figure('Color', 'w', 'Position', [179 185 1098 496]);
set(gcf, 'OuterPosition', [179 185 1098 596]);
set(gcf, 'InnerPosition', [179 185 1098 496+100]);

%--------------------------------------------------------------------------
% Case1: (PASL) GKM vs Numeric
%--------------------------------------------------------------------------
subplot(2,4,1); hold on; grid on; 
set(gca, 'Box', 'On', 'Color', background_color, 'FontSize', FontSize);
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
plot(t1, PASL_analytic1*1e2, '-', t1, PASL_numeric1*1e2, '--', 'LineWidth', LineWidth);
legend('GKM (PASL)', 'Numeric', 'Location', 'southwest');
legend('boxoff');
xlabel('Time (sec)', 'FontSize', FontSize);
title({'GKM (PASL) vs Numeric', sprintf('\\tau = %1.0f msec', dt1*1e3)}, 'FontSize', FontSize);
ax1 = gca;
hYLabel1 = ylabel('Signal / M_0 (%)', 'FontSize', FontSize);

%--------------------------------------------------------------------------
% Case2: (PASL) GKM vs Numeric
%--------------------------------------------------------------------------
subplot(2,4,2); hold on; grid on;
set(gca, 'Box', 'On', 'Color', background_color, 'FontSize', FontSize);
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
plot(t2, PASL_analytic2*1e2, '-', t2, PASL_numeric2*1e2, '--', 'LineWidth', LineWidth);
legend('GKM (PASL)', 'Numeric', 'Location', 'southwest');
legend('boxoff');
xlabel('Time (sec)', 'FontSize', FontSize);
title({'GKM (PASL) vs Numeric', sprintf('\\tau = %1.0f msec', dt2*1e3)}, 'FontSize', FontSize);
ax2 = gca;
position = get(ax2, 'Position');
set(gca, 'Position', [position(1)-spacing position(2) position(3) position(4)]);

%--------------------------------------------------------------------------
% Case3: (PCASL) GKM vs Numeric
%--------------------------------------------------------------------------
subplot(2,4,3); hold on; grid on;
set(gca, 'Box', 'On', 'Color', background_color, 'FontSize', FontSize);
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
plot(t1, PCASL_analytic1*1e2, '-', t1, PCASL_numeric1*1e2, '--', 'LineWidth', LineWidth);
legend('GKM (PCASL)', 'Numeric', 'Location', 'southwest');
legend('boxoff');
xlabel('Time (sec)', 'FontSize', FontSize);
title({'GKM (PCASL) vs Numeric', sprintf('\\tau = %1.0f msec', dt1*1e3)}, 'FontSize', FontSize);
ax3 = gca;
position = get(ax3, 'Position');
set(ax3, 'Position', [position(1)-spacing*2 position(2) position(3) position(4)]);

%--------------------------------------------------------------------------
% Case4: (PCASL) GKM vs Numeric
%--------------------------------------------------------------------------
subplot(2,4,4); hold on; grid on; 
set(gca, 'Box', 'On', 'Color', background_color, 'FontSize', FontSize);
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
plot(t2, PCASL_analytic2*1e2, '-', t2, PCASL_numeric2*1e2, '--', 'LineWidth', LineWidth);
legend('GKM (PCASL)', 'Numeric', 'Location', 'southwest');
legend('boxoff');
xlabel('Time (sec)', 'FontSize', FontSize);
title({'GKM (PCASL) vs Numeric', sprintf('\\tau = %1.0f msec', dt2*1e3)}, 'FontSize', FontSize);
ax4 = gca;
position = get(ax4, 'Position');
set(ax4, 'Position', [position(1)-spacing*3 position(2) position(3) position(4)]);

%--------------------------------------------------------------------------
% Case1: signal difference
%--------------------------------------------------------------------------
subplot(2,4,5); hold on; grid on; 
set(gca, 'Box', 'On', 'Color', background_color, 'FontSize', FontSize);
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
plot(t1, (PASL_analytic1-PASL_numeric1)*1e2, 'LineWidth', LineWidth);
xlabel('Time (sec)', 'FontSize', FontSize);
title({'Signal Difference', '(GKM - Numeric)'});
ax5 = gca;
position = get(ax5, 'Position');
set(ax5, 'Position', [position(1) position(2)-vspacing position(3) position(4)]);
hYLabel2 = ylabel('Signal / M_0 (%)', 'FontSize', FontSize);

%--------------------------------------------------------------------------
% Case2: signal difference
%--------------------------------------------------------------------------
subplot(2,4,6); hold on; grid on; 
set(gca, 'Box', 'On', 'Color', background_color, 'FontSize', FontSize);
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
plot(t2, (PASL_analytic2-PASL_numeric2)*1e2, 'LineWidth', LineWidth);
xlabel('Time (sec)', 'FontSize', FontSize);
title({'Signal Difference', '(GKM - Numeric)'});
ax6 = gca;
position = get(ax6, 'Position');
set(ax6, 'Position', [position(1)-spacing position(2)-vspacing position(3) position(4)]);

%--------------------------------------------------------------------------
% Case3: signal difference
%--------------------------------------------------------------------------
subplot(2,4,7); hold on; grid on; 
set(gca, 'Box', 'On', 'Color', background_color, 'FontSize', FontSize);
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
plot(t1, (PCASL_analytic1-PCASL_numeric1)*1e2, 'LineWidth', LineWidth);
xlabel('Time (sec)', 'FontSize', FontSize);
title({'Signal Difference', '(GKM - Numeric)'});
ax7 = gca;
position = get(gca, 'Position');
set(gca, 'Position', [position(1)-spacing*2 position(2)-vspacing position(3) position(4)]);

%--------------------------------------------------------------------------
% Case4: signal difference
%--------------------------------------------------------------------------
subplot(2,4,8); hold on; grid on;
set(gca, 'Box', 'On', 'Color', background_color, 'FontSize', FontSize);
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
plot(t2, (PCASL_analytic2-PCASL_numeric2)*1e2, 'LineWidth', LineWidth);
xlabel('Time (sec)', 'FontSize', FontSize);
title({'Signal Difference', '(GKM - Numeric)'});
ax8 = gca;
position = get(ax8, 'Position');
set(ax8, 'Position', [position(1)-spacing*3 position(2)-vspacing position(3) position(4)]);

% Set ylim for all subplots
y_max = max(abs(get(ax4, 'YLim')));
ylimits = [0 y_max]*1.1;
set(ax1, 'YLim', ylimits);
set(ax2, 'YLim', ylimits);
set(ax3, 'YLim', ylimits);
set(ax4, 'YLim', ylimits);

% Set ylim for all subplots
y_max = max(abs(get(ax6, 'YLim')));
ylimits = [-y_max y_max];
set(ax5, 'YLim', ylimits, 'YTick', (-y_max:0.02:y_max), 'YTickLabel', (-y_max:0.02:y_max));
set(ax6, 'YLim', ylimits, 'YTick', (-y_max:0.02:y_max), 'YTickLabel', (-y_max:0.02:y_max));
set(ax7, 'YLim', ylimits, 'YTick', (-y_max:0.02:y_max), 'YTickLabel', (-y_max:0.02:y_max));
set(ax8, 'YLim', ylimits, 'YTick', (-y_max:0.02:y_max), 'YTickLabel', (-y_max:0.02:y_max));

%--------------------------------------------------------------------------
% Labeling
%--------------------------------------------------------------------------
labeling_width      = 0.2;
labeling_color      = [204 204 204] / 255;
labeling_text_color = [128 128 128] / 255;

YLim = get(ax1, 'YLim');
x = [t_labeling t_labeling t_labeling+labeling_width t_labeling+labeling_width];
y = [0 YLim(2) YLim(2) 0];

h1 = patch(ax1, 'XData', x, 'YData', y, 'FaceColor', labeling_color, 'EdgeColor', 'none', 'FaceAlpha', 1, 'LineWidth', 1, 'LineStyle', '-');
set(h1, 'HandleVisibility', 'Off');
text(ax1, t_labeling+0.2, YLim(2), 'RF labeling', 'Color', labeling_text_color, 'VerticalAlignment', 'top');

h2 = patch(ax2, 'XData', x, 'YData', y, 'FaceColor', labeling_color, 'EdgeColor', 'none', 'FaceAlpha', 1, 'LineWidth', 1, 'LineStyle', '-');
set(h2, 'HandleVisibility', 'Off');
text(ax2, t_labeling+0.2, YLim(2), 'RF labeling', 'Color', labeling_text_color, 'VerticalAlignment', 'top');

h3 = patch(ax3, 'XData', x, 'YData', y, 'FaceColor', labeling_color, 'EdgeColor', 'none', 'FaceAlpha', 1, 'LineWidth', 1, 'LineStyle', '-');
set(h3, 'HandleVisibility', 'Off');
text(ax3, t_labeling+0.2, YLim(2), 'RF labeling', 'Color', labeling_text_color, 'VerticalAlignment', 'top');

h4 = patch(ax4, 'XData', x, 'YData', y, 'FaceColor', labeling_color, 'EdgeColor', 'none', 'FaceAlpha', 1, 'LineWidth', 1, 'LineStyle', '-');
set(h4, 'HandleVisibility', 'Off');
text(ax4, t_labeling+0.2, YLim(2), 'RF labeling', 'Color', labeling_text_color, 'VerticalAlignment', 'top');

text(-0.8-5.15*3,  0.108, '(a)', 'FontWeight', 'bold', 'FontSize', 12);
text(-0.8-5.15*2,  0.108, '(b)', 'FontWeight', 'bold', 'FontSize', 12);
text(-5.95      ,  0.108, '(c)', 'FontWeight', 'bold', 'FontSize', 12);
text(-0.8       ,  0.108, '(d)', 'FontWeight', 'bold', 'FontSize', 12);

text(-0.8-5.15*3, -0.085, '(e)', 'FontWeight', 'bold', 'FontSize', 12);
text(-0.8-5.15*2, -0.085, '(f)', 'FontWeight', 'bold', 'FontSize', 12);
text(-5.95      , -0.085, '(g)', 'FontWeight', 'bold', 'FontSize', 12);
text(-0.8       , -0.085, '(h)', 'FontWeight', 'bold', 'FontSize', 12);

%% Adjust the location of the ylabel in the first row
pause(1);
position1 = get(hYLabel1, 'Position');
position2 = get(hYLabel2, 'Position');
set(hYLabel1, 'Position', [position2(1) position1(2) position1(3)]);

%% Save as a .tiff file
export_fig('Figure2', '-r864', '-tif');
