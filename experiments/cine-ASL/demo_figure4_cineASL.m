% demo_figure4_cineASL.m
% Written by Namgyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 10/05/2019, Last modified: 10/06/2019

%% Clean slate
close all; clear all; clc;

%% Add paths
addpath(genpath('D:\ASL_project\mfiles_nam\Bloch_Flow_MT'));

%% Define simulation parameters
T1      = 1.4;          % spin-lattice relaxation time of myocardium at 3T [sec]
T2      = 25e-3;        % transverse relaxation time [sec]
alpha   = 8 * pi / 180; % flip angle [rad]
F       = 6;            % rat brain blood flow [mL/g/min]
lambda  = 0.95;         % brain/blood partition coefficient [mL of blood /g of brain tissue]
beta    = 0.5;          % effective labeling efficiency of the experiment
M0      = 1;            % equilibrium z magnetization [magnetization/g tissue]
M0b     = M0 / lambda;  % M0 arterial blood = M0 tissue / lambda [magnetization/mL of blood]
RR      = 130e-3;       % interval between heart beats
TR      = 10e-3;        % repetition time [sec] 
Nlines  = 6;            % number of k-space lines
Nechoes = 12;           % number of echoes
alpha0  = 2 * beta;     % labelling efficiency; 1 for saturation, 2 for inversion, and 0 for control

%% Compare cine-ASL to Bloch-Flow
t_blochFlow_all = cell(4,1);

mx_blochFlow_all = cell(4,1);
my_blochFlow_all = cell(4,1);
mz_blochFlow_all = cell(4,1);

t_A_t_all = cell(4,1);
t_B_t_all = cell(4,1);
t_A_c_all = cell(4,1);
t_B_c_all = cell(4,1);

Mt_A_all = cell(4,1);
Mt_B_all = cell(4,1);
Mc_A_all = cell(4,1);
Mc_B_all = cell(4,1);

RDs = [0.01; 0.01; 3.5 ; 3.5 ]; % recovery delay [sec]
tps = [0.26; 3.12; 0.26; 3.12]; % the duration of the acquisition phase [sec]

for idx1 = 1:4
    RD = RDs(idx1); % recovery delay [sec]
    tp = tps(idx1); % duration of the acquisition phase [sec]: tp = Ncine * (Nechoes + 1) * TR

    %% Calculate the number of samples
    N_A_t = ceil(tp / TR);
    N_B_t = ceil(RD / TR);
    N_A_c = ceil(tp / TR);
    N_B_c = ceil(RD / TR);
    N     = Nlines * (N_A_t + N_B_t + N_A_c + N_B_c);

    %% Calculate b1, gradient and taus
    gam = 4257.746778 * 2 * pi * 1e-2; % [rad/sec/uT]

    b1 = zeros(N,1, 'double'); % [uT]
    gr = zeros(N,1, 'double'); % [G/cm]

    phaseA_tag_index  = reshape(bsxfun(@plus,                         (1:N_A_t).', (N_A_t + N_B_t + N_A_c + N_B_c) * (0:Nlines-1)), [N_A_t*Nlines 1]);
    phaseB_tag_index  = reshape(bsxfun(@plus, N_A_t                 + (1:N_B_t).', (N_A_t + N_B_t + N_A_c + N_B_c) * (0:Nlines-1)), [N_B_t*Nlines 1]);
    phaseA_ctrl_index = reshape(bsxfun(@plus, N_A_t + N_B_t         + (1:N_A_c).', (N_A_t + N_B_t + N_A_c + N_B_c) * (0:Nlines-1)), [N_A_c*Nlines 1]);
    phaseB_ctrl_index = reshape(bsxfun(@plus, N_A_t + N_B_t + N_A_c + (1:N_B_c).', (N_A_t + N_B_t + N_A_c + N_B_c) * (0:Nlines-1)), [N_B_c*Nlines 1]);

    phase_incr = pi; % phase cycling increment [rad]
    rf_phase_cycling = (0:N_A_t-1).' * phase_incr; % [rad]

    % [rad] / [sec] / ([rad/sec/uT]) => [uT]
    b1(phaseA_tag_index)  = repmat(alpha / TR / gam * exp(1j * rf_phase_cycling), [Nlines 1]); % [uT]
    b1(phaseA_ctrl_index) = repmat(alpha / TR / gam * exp(1j * rf_phase_cycling), [Nlines 1]); % [uT]

    taus = TR * ones(N,1, 'double');  % [sec]
    t_blochFlow = (0:N-1).' * TR; % [sec]

    %% Calculate ASL_bolus [magnetization/g tissue/sec]
    ASL_bolus = zeros(N,1, 'double');
    ASL_bolus(phaseA_tag_index) = -(F / lambda / 60) * M0 * alpha0;

    %% Perform Bloch-Flow simulation
    % Get the magnetization immediately after an RF pulse
    tic; [mx_blochFlow,my_blochFlow,mz_blochFlow] = bloch_flow_spoiling(b1*1e-2, gr, taus, 0*taus, T1, T2, 0, 0, 2, 0, 0, 1, F, lambda, M0, ASL_bolus); toc;

    t_blochFlow_all{idx1}  = t_blochFlow;
    mx_blochFlow_all{idx1} = mx_blochFlow;
    my_blochFlow_all{idx1} = my_blochFlow;
    mz_blochFlow_all{idx1} = mz_blochFlow;

    %% Calculate the steady-state magnetization for cine-FLASH
    Mss = M0 * (1 - exp(-TR / T1)) / (1 - cos(alpha) * exp(-TR / T1));

    %% Calculate the decay time constants
    T1_star    = 1 / (1 / T1 - log(cos(alpha)) / TR);                     % [sec]
    T1app_star = 1 / (1 / T1 + (F / lambda / 60) - log(cos(alpha)) / TR); % [sec], during phase A
    T1app      = 1 / (1 / T1 + (F / lambda / 60));                        % [sec], used during phase B

    %% Calculate Ma_ctrl and Ma_tag
    Ma_ctrl = M0b;                  % [magnetization/mL of blood]
    Ma_tag  = (1 - 2 * beta) * M0b; % [magnetization/mL of blood]

    %% Calculate the asymptotic value of the magnetization at infinite time
    % [mL/g/min] * [magnetization/mL of blood] * [min/60sec] => [magnetization/g tissue/sec]
    Minf_A_t = Mss * T1app_star / T1_star + (F * Ma_tag  / 60) * T1app_star;
    Minf_A_c = Mss * T1app_star / T1_star + (F * Ma_ctrl / 60) * T1app_star;

    %% Calculate the theoretical signal evolution for cine-ASL
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

        %------------------------------------------------------------------
        % Phase A, tag
        %------------------------------------------------------------------
        t_A_t(:,n+1) = (0:nt_A-1).' * TR + 2 * n * (tp + RD); % [sec]
        Mt_A(:,n+1)  = Minf_A_t + (MA_t - Minf_A_t) * exp(-(t_A_t(:,n+1) - 2 * n * (tp + RD)) / T1app_star);

        %------------------------------------------------------------------
        % Phase B, tag
        %------------------------------------------------------------------
        t_B_t(:,n+1) = (0:nt_B-1).' * TR + 2 * n * (tp + RD) + tp; % [sec]
        Mt_B(:,n+1)  = M0 + (MB_t - M0) * exp(-(t_B_t(:,n+1) - (2 * n * (tp + RD) + tp)) / T1app);

        %------------------------------------------------------------------
        % Phase A, control
        %------------------------------------------------------------------
        t_A_c(:,n+1) = (0:nt_A-1).' * TR + (2 * n + 1) * (tp + RD); % [sec]
        Mc_A(:,n+1)  = Minf_A_c + (MA_c - Minf_A_c) * exp(-(t_A_c(:,n+1) - (2 * n + 1) * (tp + RD)) / T1app_star);

        %------------------------------------------------------------------
        % Phase B, control
        %------------------------------------------------------------------
        t_B_c(:,n+1) = (0:nt_B-1).' * TR + (2 * n + 1) * (tp + RD) + tp; % [sec]
        Mc_B(:,n+1)  = M0 + (MB_c - M0) * exp(-(t_B_c(:,n+1) - ((2 * n + 1) * (tp + RD) + tp)) / T1app);

        %------------------------------------------------------------------
        % Update MA_t and MB_t
        %------------------------------------------------------------------
        MA_t = yA_t + x * MA_c;
        MB_t = yB_t + x * MB_c;
    end

    t_A_t_all{idx1} = t_A_t;
    t_B_t_all{idx1} = t_B_t;
    t_A_c_all{idx1} = t_A_c;
    t_B_c_all{idx1} = t_B_c;

    Mt_A_all{idx1} = Mt_A;
    Mt_B_all{idx1} = Mt_B;
    Mc_A_all{idx1} = Mc_A;
    Mc_B_all{idx1} = Mc_B;
end

%% Display results
figure('Color', 'w', 'Position', [2 124 894 692]);
color_order = get(gca, 'colororder');
vec = @ (x) x(:);

%--------------------------------------------------------------------------
% Figure 1 (a)
%--------------------------------------------------------------------------
subplot(2,2,1); hold on; set(gca, 'Box', 'On'); grid on;
t_cineASL  = vec(cat(1, t_A_t_all{1}, t_B_t_all{1}, t_A_c_all{1}, t_B_c_all{1}));
mz_cineASL = vec(cat(1, Mt_A_all{1} , Mt_B_all{1} , Mc_A_all{1} , Mc_B_all{1}));
t_blochFlow = t_blochFlow_all{1};
mz_blochFlow = mz_blochFlow_all{1};
plot(t_cineASL  , mz_cineASL  , 'Color', color_order(1,:), 'LineWidth', 1);
plot(t_blochFlow, mz_blochFlow, 'Color', color_order(2,:), 'LineWidth', 1);
xlim([0 3.5]);
ylim([0.3 1]);
legend_text = {'Cine-ASL', 'Numeric'};
legend(legend_text);
xlabel('Time (sec)');
ylabel('Signal / M_0');
title(sprintf('Nlines = %d, tp = %4.2f sec, RD = %4.2f sec', Nlines, tps(1), RDs(1)));
difference1 = mz_cineASL - mz_blochFlow;

%--------------------------------------------------------------------------
% Figure 1 (b)
%--------------------------------------------------------------------------
subplot(2,2,2); hold on; set(gca, 'Box', 'On'); grid on;
t_cineASL  = vec(cat(1, t_A_t_all{2}, t_B_t_all{2}, t_A_c_all{2}, t_B_c_all{2}));
mz_cineASL = vec(cat(1, Mt_A_all{2} , Mt_B_all{2} , Mc_A_all{2} , Mc_B_all{2}));
t_blochFlow = t_blochFlow_all{2};
mz_blochFlow = mz_blochFlow_all{2};
plot(t_cineASL  , mz_cineASL  , 'Color', color_order(1,:), 'LineWidth', 1);
plot(t_blochFlow, mz_blochFlow, 'Color', color_order(2,:), 'LineWidth', 1);
xlim([0 40]);
ylim([0.3 1]);
legend(legend_text);
xlabel('Time (sec)');
ylabel('Signal / M_0');
title(sprintf('Nlines = %d, tp = %4.2f sec, RD = %4.2f sec', Nlines, tps(2), RDs(2)));
difference2 = mz_cineASL - mz_blochFlow;

%--------------------------------------------------------------------------
% Figure 1 (c)
%--------------------------------------------------------------------------
subplot(2,2,3); hold on; set(gca, 'Box', 'On'); grid on;
t_cineASL  = vec(cat(1, t_A_t_all{3}, t_B_t_all{3}, t_A_c_all{3}, t_B_c_all{3}));
mz_cineASL = vec(cat(1, Mt_A_all{3} , Mt_B_all{3} , Mc_A_all{3} , Mc_B_all{3}));
t_blochFlow = t_blochFlow_all{3};
mz_blochFlow = mz_blochFlow_all{3};
plot(t_cineASL  , mz_cineASL  , 'Color', color_order(1,:), 'LineWidth', 1);
plot(t_blochFlow, mz_blochFlow, 'Color', color_order(2,:), 'LineWidth', 1);
xlim([0 50]);
ylim([0.3 1]);
legend(legend_text, 'Location', 'Southeast');
xlabel('Time (sec)');
ylabel('Signal / M_0');
title(sprintf('Nlines = %d, tp = %4.2f sec, RD = %4.2f sec', Nlines, tps(3), RDs(3)));
difference3 = mz_cineASL - mz_blochFlow;

%--------------------------------------------------------------------------
% Figure 1 (d)
%--------------------------------------------------------------------------
subplot(2,2,4); hold on; set(gca, 'Box', 'On'); grid on;
t_cineASL  = vec(cat(1, t_A_t_all{4}, t_B_t_all{4}, t_A_c_all{4}, t_B_c_all{4}));
mz_cineASL = vec(cat(1, Mt_A_all{4} , Mt_B_all{4} , Mc_A_all{4} , Mc_B_all{4}));
t_blochFlow = t_blochFlow_all{4};
mz_blochFlow = mz_blochFlow_all{4};
plot(t_cineASL  , mz_cineASL  , 'Color', color_order(1,:), 'LineWidth', 1);
plot(t_blochFlow, mz_blochFlow, 'Color', color_order(2,:), 'LineWidth', 1);
xlim([0 80]);
ylim([0.3 1]);
legend(legend_text, 'Location', 'Southeast');
xlabel('Time (sec)');
ylabel('Signal / M_0');
title(sprintf('Nlines = %d, tp = %4.2f sec, RD = %4.2f sec', Nlines, tps(4), RDs(4)));
difference4 = mz_cineASL - mz_blochFlow;

text(-12-105, 0.21+0.98, '(a)', 'FontWeight', 'bold', 'FontSize', 12);
text(-12    , 0.21+0.98, '(b)', 'FontWeight', 'bold', 'FontSize', 12);
text(-12-105, 0.21     , '(c)', 'FontWeight', 'bold', 'FontSize', 12);
text(-12    , 0.21     , '(d)', 'FontWeight', 'bold', 'FontSize', 12);

figure;
subplot(2,2,1);
plot(difference1);

subplot(2,2,2);
plot(difference2);

subplot(2,2,3);
plot(difference3);

subplot(2,2,4);
plot(difference4);


%% Save as a .tiff file
%export_fig('comparison_cineASL_vs_Bloch_flow', '-r864', '-tif');
