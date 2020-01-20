% demo_Capron_2013_MRM.m
% Written by Namgyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 10/04/2019, Last modified: 10/05/2019

%% Clean slate
close all; clear all; clc;

%% Add paths
addpath(genpath('D:\ASL_project\mfiles_nam\Bloch-Flow'));

%% Define simulation parameters
T1     = 1.4;          % spin-lattice relaxation time of myocardium at 3T [sec]
alpha  = 8 * pi / 180; % flip angle [rad]
F      = 6;            % rat brain blood flow [mL/g/min]
lambda = 0.95;         % brain/blood partition coefficient [mL of blood /g of brain tissue]
beta   = 0.5;          % effective labeling efficiency of the experiment
M0     = 1;            % equilibrium z magnetization [magnetization/g tissue]
M0b    = M0 / lambda;  % M0 arterial blood = M0 tissue / lambda [magnetization/mL of blood]
RR     = 130e-3;       % interval between heart beats
TR     = 10e-3;        % repetition time [sec] 

%% Calculate the steady-state magnetization for cine-FLASH
Mss = M0 * (1 - exp(-TR / T1)) / (1 - cos(alpha) * exp(-TR / T1));

%% Calculate the decay time constants
T1_star    = 1 / (1 / T1 - log(cos(alpha)) / TR); % [sec]
T1app_star = 1 / (1 / T1 + (F / lambda / 60) - log(cos(alpha)) / TR); % [sec], during phase A 
T1app      = 1 / (1 / T1 + (F / lambda / 60)); % [sec], used during phase B 

%% Calculate Ma_ctrl and Ma_tag
Ma_ctrl = M0b;                  % [magnetization/mL of blood]
Ma_tag  = (1 - 2 * beta) * M0b; % [magnetization/mL of blood]

%% Calculate the asymptotic value of the magnetization at infinite time
% [mL/g/min] * [magnetization/mL of blood] * [min/60sec] => [magnetization/g tissue/sec]
Minf_A_t = Mss * T1app_star / T1_star + (F * Ma_tag  / 60) * T1app_star;
Minf_A_c = Mss * T1app_star / T1_star + (F * Ma_ctrl / 60) * T1app_star;

%%
Nlines  = 6;    % number of k-space lines
Ncine   = 24;   % number of cine blocks
Nechoes = 12;   % number of echoes
tp      = 3.12; % duration of the acquisition phase [sec]: tp = Ncine * (Nechoes + 1) * TR
RD      = 0.01; % recovery delay [sec]


yB_t = M0 * (1 - exp(-RD / T1app)) * exp(-tp / T1app_star) + Minf_A_t * (1 - exp(-tp / T1app_star));
yB_c = M0 * (1 - exp(-RD / T1app)) * exp(-tp / T1app_star) + Minf_A_c * (1 - exp(-tp / T1app_star));
yA_t = M0 * (1 - exp(-RD / T1app)) + Minf_A_c * (1 - exp(-tp / T1app_star)) * exp(-RD / T1app);
yA_c = M0 * (1 - exp(-RD / T1app)) + Minf_A_t * (1 - exp(-tp / T1app_star)) * exp(-RD / T1app);

x = exp(-RD / T1app) * exp(-tp / T1app_star);


nt_A = ceil(tp / TR);
nt_B = ceil(RD / TR);

% Initialize MA_t and MB_t
MA_t = M0;
MB_t = Minf_A_t + (MA_t - Minf_A_t) * exp(-tp / T1app_star); % [A1]

t_A_t = zeros(nt_A,Nlines, 'double');
t_B_t = zeros(nt_B,Nlines, 'double');
t_A_c = zeros(nt_A,Nlines, 'double');
t_B_c = zeros(nt_B,Nlines, 'double');

Mt_A = zeros(nt_A,Nlines, 'double');
Mt_B = zeros(nt_B,Nlines, 'double');
Mc_A = zeros(nt_A,Nlines, 'double');
Mc_B = zeros(nt_B,Nlines, 'double');

for n = 0:Nlines-1
    MA_c = yA_c + x * MA_t;
    MB_c = yB_c + x * MB_t;

    %----------------------------------------------------------------------
    % Phase A, tag
    %----------------------------------------------------------------------
    t_A_t(:,n+1) = (0:nt_A-1).' * TR + 2 * n * (tp + RD); % [sec]
    Mt_A(:,n+1)  = Minf_A_t + (MA_t - Minf_A_t) * exp(-(t_A_t(:,n+1) - 2 * n * (tp + RD)) / T1app_star);

    %----------------------------------------------------------------------
    % Phase B, tag
    %----------------------------------------------------------------------
    t_B_t(:,n+1) = (0:nt_B-1).' * TR + 2 * n * (tp + RD) + tp; % [sec]
    Mt_B(:,n+1)  = M0 + (MB_t - M0) * exp(-(t_B_t(:,n+1) - (2 * n * (tp + RD) + tp)) / T1app);

    %----------------------------------------------------------------------
    % Phase A, control
    %----------------------------------------------------------------------
    t_A_c(:,n+1) = (0:nt_A-1).' * TR + (2 * n + 1) * (tp + RD); % [sec]
    Mc_A(:,n+1)  = Minf_A_c + (MA_c - Minf_A_c) * exp(-(t_A_c(:,n+1) - (2 * n + 1) * (tp + RD)) / T1app_star);

    %----------------------------------------------------------------------
    % Phase B, control
    %----------------------------------------------------------------------
    t_B_c(:,n+1) = (0:nt_B-1).' * TR + (2 * n + 1) * (tp + RD) + tp; % [sec]
    Mc_B(:,n+1)  = M0 + (MB_c - M0) * exp(-(t_B_c(:,n+1) - ((2 * n + 1) * (tp + RD) + tp)) / T1app);

    %----------------------------------------------------------------------
    % Update MA_t and MB_t
    %----------------------------------------------------------------------
    MA_t = yA_t + x * MA_c;
    MB_t = yB_t + x * MB_c;
end

%%
figure('Color', 'w');
color_order = get(gca, 'colororder');

hold on; set(gca, 'Box', 'On'); grid on;
for n = 0:Nlines-1
    plot(t_A_t(:,n+1), Mt_A(:,n+1), '-' , 'Color', color_order(1,:), 'LineWidth', 1);
    plot(t_B_t(:,n+1), Mt_B(:,n+1), '.-', 'Color', color_order(2,:), 'LineWidth', 1, 'MarkerSize', 10);
    plot(t_A_c(:,n+1), Mc_A(:,n+1), '-' , 'Color', color_order(3,:), 'LineWidth', 1);
    plot(t_B_c(:,n+1), Mc_B(:,n+1), '.-', 'Color', color_order(4,:), 'LineWidth', 1, 'MarkerSize', 10);
end
ylim([0.3 1]);
legend('Phase A, tag', 'Phase B, tag', 'Phase A, ctrl', 'Phase B, ctrl');
xlabel('Time (sec)');
ylabel('Signal / M_0');
title(sprintf('Nlines = %d, tp = %4.2f sec, RD = %4.2f sec', Nlines, tp, RD));

% %% Save as a .tiff file
% %export_fig(sprintf('Figure2_%s', result_text), '-r864', '-tiff');
% end
