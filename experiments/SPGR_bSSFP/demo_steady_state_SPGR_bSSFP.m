% demo_steady_state_SPGR_bSSFP.m
% Written by Namgyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 11/18/2019, Last modified: 11/23/2019

%% Clean slate
close all; clear all; clc;

%% Add paths
addpath D:\ASL_project\mfiles_nam\Bloch-Flow\lib;
addpath D:\ASL_project\mfiles_nam\Bloch-Flow\experiments\SPGR;
addpath D:\ASL_project\mfiles_nam\Bloch-Flow\experiments\Gloor_2008_MRM;

%% Define parameters for white matter at 1.5T
M0f    = 1;          % equilibrium magnetiztion of the free pool
R1f    = 1 / 585e-3; % Free pool longitudinal relaxation rate    [1/sec]
R2f    = 1 / 81e-3;  % Free pool transverse relaxation rate      [1/sec]
f      = 0.157;      % the fractional size of the semisolid pool
M0s    = f * M0f;    % equilibrium magnetiztion of the semisolid pool
kf     = 4.45;       % pseudo-first-order rate constant          [1/sec]
k      = kf / M0s;   % fundamental rate constant                 [1/sec]
R1s    = 1;          % Semisolid longitudinal relaxation rate    [1/sec]
T2s    = 12e-6;      % Semisolid transverse relaxation time      [sec]
G      = 14e-6;      % lineshape value at zero frequency G(0)    [sec]
F      = 4;          % brain blood flow [mL/g/min] WM: 34 [ml/100 ml/min], GM: 42 [ml/100 ml/min]
lambda = 0.9;        % brain/blood partition coefficient [mL of blood /g of brain tissue]

dt        = 4e-6;              % RF dwell time [sec]
tau_bssfp = 230e-6;            % pulse duration [sec] (T_RF in the paper)
TR_bssfp  = 2.92e-3;           % repetition time [sec]
TE_bssfp  = TR_bssfp / 2;      % echo time [sec]
tau_spgr  = 2e-3;              % pulse duration [sec] (T_RF in the paper)
TR_spgr   = 5e-3;              % repetition time [sec]
TE_spgr   = TR_spgr / 2;       % echo time [sec]
flips = (1:1:80).' * pi / 180; % flip angles [rad]


% %%
% % % Malik's numerical simulation
% M0f    = 0.853;      % equilibrium magnetiztion of the free pool
% R1f    = 1.54;       % Free pool longitudinal relaxation rate    [1/sec]
% R2f    = 12.5;       % Free pool transverse relaxation rate      [1/sec]
% M0s    = 0.147;      % equilibrium magnetiztion of the semisolid pool
% f      = M0s / M0f;  % the fractional size of the semisolid pool
% k      = 65;         % fundamental rate constant                 [1/sec]
% R1s    = 1;          % Semisolid longitudinal relaxation rate    [1/sec]
% T2s    = 12.5e-6;    % Semisolid transverse relaxation time      [sec]
% G      = 14e-6;      % lineshape value at zero frequency G(0)    [sec]
% F      = 4;          % brain blood flow [mL/g/min] WM: 34 [ml/100 ml/min], GM: 42 [ml/100 ml/min]
% lambda = 0.9;        % brain/blood partition coefficient [mL of blood /g of brain tissue]
% 
% dt = 4e-6;           % RF dwell time [sec]
% tau = 2e-3;          % pulse duration [sec] (T_RF in the paper)
% TR  = 5e-3;          % repetition time [sec]
% TE  = TR / 2;        % echo time [sec]
% flips = (1:1:80).' * pi / 180; % flip angles [rad]

tissuepars = struct;
tissuepars.free.M0 = M0f;
tissuepars.free.R1 = R1f;
tissuepars.free.R2 = R2f;
tissuepars.semi.M0 = M0s;
tissuepars.semi.R1 = R1s;
tissuepars.semi.T2 = T2s;
tissuepars.k       = k;
tissuepars.f       = f;
tissuepars.G       = G;

tissuepars_noMT = struct;
tissuepars_noMT.free.M0 = M0f;
tissuepars_noMT.free.R1 = R1f;
tissuepars_noMT.free.R2 = R2f;
tissuepars_noMT.semi.M0 = 0;
tissuepars_noMT.semi.R1 = R1s;
tissuepars_noMT.semi.T2 = T2s;
tissuepars_noMT.k       = 0;
tissuepars_noMT.f       = 0;
tissuepars_noMT.G       = G;

%% Calculate the time axis
nt_bssfp = ceil(tau_bssfp / dt); % number of samples in time domain
t_bssfp = (0:nt_bssfp-1).' * dt; % [sec]

nt_spgr = ceil(tau_spgr / dt); % number of samples in time domain
t_spgr = (0:nt_spgr-1).' * dt; % [sec]

%% Calculate a Hamming windowed sinc pulse
TBW = 2.2; % time-bandwidth product [sec] * [Hz]
pulse_bssfp = wsinc(TBW, nt_bssfp);
pulse_bssfp = pulse_bssfp / max(pulse_bssfp);

pulse_spgr = wsinc(TBW, nt_spgr);
pulse_spgr = pulse_spgr / max(pulse_spgr);

%% Create RF pulses
% 1T = 1e4G = 1e6uT => 1G = 1e2uT
% Gamma for current nucleus [Hz/G] * [2*pi rad/cycle] * [G/1e2 uT] => [rad/sec/uT]
gam = 4257.59 * 2 * pi * 1e-2; % [rad/sec/uT]
nr_flips = length(flips);
b1sqrd_bssfp = zeros(nr_flips,1, 'double'); % [uT^2]
b1sqrd_spgr  = zeros(nr_flips,1, 'double'); % [uT^2]

for idx = 1:nr_flips
    %----------------------------------------------------------------------
    % Calculate the amplitude of an RF pulse [uT]
    % Note: The area of an RF pulse is a flip angle
    % flipangle = gam * b1 * sum(pulse) * dt
    % b1 = flipangle / (gam * sum(pulse) * dt)
    % [rad] / ([rad/sec/uT] * [sec]) => [uT]
    %----------------------------------------------------------------------
    b1_max = flips(idx) / (gam * sum(pulse_bssfp) * dt); % [uT]
    pulseSB_bssfp = pulse_bssfp * b1_max; % [uT]
    b1sqrd_bssfp(idx,1) = sum(abs(pulseSB_bssfp).^2) * dt / tau_bssfp; % [uT^2]

    b1_max = flips(idx) / (gam * sum(pulse_spgr) * dt); % [uT]
    pulseSB_spgr = pulse_spgr * b1_max; % [uT]
    b1sqrd_spgr(idx,1) = sum(abs(pulseSB_spgr).^2) * dt / tau_spgr; % [uT^2]
end

%% Calculate the mean saturation rate [rad/sec]
% [rad/sec/uT]^2 * [uT^2] * [sec] => [rad/sec]
W_bssfp = pi * gam^2 * b1sqrd_bssfp * G; % [rad/sec]
W_spgr  = pi * gam^2 * b1sqrd_spgr  * G; % [rad/sec]

%% Calculate the steady-state magnetization
bssfp_approx_noMT_noflow = zeros(4,nr_flips, 'double');
bssfp_approx_noMT_flow   = zeros(4,nr_flips, 'double');
bssfp_approx_MT_noflow   = zeros(4,nr_flips, 'double');
bssfp_approx_MT_flow     = zeros(4,nr_flips, 'double');

bssfp_exact_noMT_noflow = zeros(4,nr_flips, 'double');
bssfp_exact_noMT_flow   = zeros(4,nr_flips, 'double');
bssfp_exact_MT_noflow   = zeros(4,nr_flips, 'double');
bssfp_exact_MT_flow     = zeros(4,nr_flips, 'double');

spgr_approx_noMT_noflow = zeros(4,nr_flips, 'double');
spgr_approx_noMT_flow   = zeros(4,nr_flips, 'double');
spgr_approx_MT_noflow   = zeros(4,nr_flips, 'double');
spgr_approx_MT_flow     = zeros(4,nr_flips, 'double');

spgr_exact_noMT_noflow = zeros(4,nr_flips, 'double');
spgr_exact_noMT_flow   = zeros(4,nr_flips, 'double');
spgr_exact_MT_noflow   = zeros(4,nr_flips, 'double');
spgr_exact_MT_flow     = zeros(4,nr_flips, 'double');

tic;
for idx = 1:nr_flips
    %----------------------------------------------------------------------
    % bSSFP
    %----------------------------------------------------------------------
    bssfp_approx_noMT_noflow(:,idx) = ssSSFP_BMF(flips(idx), 0                , TR_bssfp, TE_bssfp, tau_bssfp, 0, pi, 0, lambda, tissuepars_noMT);
    bssfp_approx_noMT_flow(:,idx)   = ssSSFP_BMF(flips(idx), 0                , TR_bssfp, TE_bssfp, tau_bssfp, 0, pi, F, lambda, tissuepars_noMT);
    bssfp_approx_MT_noflow(:,idx)   = ssSSFP_BMF(flips(idx), b1sqrd_bssfp(idx), TR_bssfp, TE_bssfp, tau_bssfp, 0, pi, 0, lambda, tissuepars);
    bssfp_approx_MT_flow(:,idx)     = ssSSFP_BMF(flips(idx), b1sqrd_bssfp(idx), TR_bssfp, TE_bssfp, tau_bssfp, 0, pi, F, lambda, tissuepars);

    bssfp_exact_noMT_noflow(:,idx)  = ssSSFP_BMF_expm(flips(idx), 0                , TR_bssfp, TE_bssfp, tau_bssfp, 0, pi, 0, lambda, tissuepars_noMT);
    bssfp_exact_noMT_flow(:,idx)    = ssSSFP_BMF_expm(flips(idx), 0                , TR_bssfp, TE_bssfp, tau_bssfp, 0, pi, F, lambda, tissuepars_noMT);
    bssfp_exact_MT_noflow(:,idx)    = ssSSFP_BMF_expm(flips(idx), b1sqrd_bssfp(idx), TR_bssfp, TE_bssfp, tau_bssfp, 0, pi, 0, lambda, tissuepars);
    bssfp_exact_MT_flow(:,idx)      = ssSSFP_BMF_expm(flips(idx), b1sqrd_bssfp(idx), TR_bssfp, TE_bssfp, tau_bssfp, 0, pi, F, lambda, tissuepars);

    %----------------------------------------------------------------------
    % SPGR
    %----------------------------------------------------------------------
    spgr_approx_noMT_noflow(:,idx) = ssSPGR_BMF(flips(idx), 0               , TR_spgr, tau_spgr, 0, 0, lambda, tissuepars_noMT);
    spgr_approx_noMT_flow(:,idx)   = ssSPGR_BMF(flips(idx), 0               , TR_spgr, tau_spgr, 0, F, lambda, tissuepars_noMT);
    spgr_approx_MT_noflow(:,idx)   = ssSPGR_BMF(flips(idx), b1sqrd_spgr(idx), TR_spgr, tau_spgr, 0, 0, lambda, tissuepars);
    spgr_approx_MT_flow(:,idx)     = ssSPGR_BMF(flips(idx), b1sqrd_spgr(idx), TR_spgr, tau_spgr, 0, F, lambda, tissuepars);

    spgr_exact_noMT_noflow(:,idx)  = ssSPGR_BMF_expm(flips(idx), 0               , TR_spgr, tau_spgr, 0, 0, lambda, tissuepars_noMT);
    spgr_exact_noMT_flow(:,idx)    = ssSPGR_BMF_expm(flips(idx), 0               , TR_spgr, tau_spgr, 0, F, lambda, tissuepars_noMT);
    spgr_exact_MT_noflow(:,idx)    = ssSPGR_BMF_expm(flips(idx), b1sqrd_spgr(idx), TR_spgr, tau_spgr, 0, 0, lambda, tissuepars);
    spgr_exact_MT_flow(:,idx)      = ssSPGR_BMF_expm(flips(idx), b1sqrd_spgr(idx), TR_spgr, tau_spgr, 0, F, lambda, tissuepars);
end
toc;

%% Calculate the single-pool bSSFP signal
bssfp_singlepool = M0f * freeman_hill(flips*180/pi, TR_bssfp, R1f, R2f);

%% Calculate the single-pool SPGR signal
spgr_singlepool = M0f * ernst(flips*180/pi, TR_spgr, R1f);

%% Display results
FontSize = 16;

figure('Color', 'w', 'Position', [7 349 1184 449]);
color_order = get(gca, 'colororder');

subplot(1,2,1); hold on; set(gca, 'Box', 'On'); grid on;
plot(flips*180/pi, bssfp_exact_noMT_noflow(2,:), '-', 'Color', color_order(1,:), 'LineWidth', 1); 
plot(flips*180/pi, bssfp_exact_noMT_flow(2,:)  , '-', 'Color', color_order(2,:), 'LineWidth', 1); 
plot(flips*180/pi, bssfp_exact_MT_noflow(2,:)  , '-', 'Color', color_order(3,:), 'LineWidth', 1); 
plot(flips*180/pi, bssfp_exact_MT_flow(2,:)    , '-', 'Color', color_order(4,:), 'LineWidth', 1); 

plot(flips(5:5:end)*180/pi, bssfp_approx_noMT_noflow(2,5:5:end), '.', 'Color', color_order(1,:), 'LineWidth', 1, 'MarkerSize', 14); 
plot(flips(5:5:end)*180/pi, bssfp_approx_noMT_flow(2,5:5:end)  , '.', 'Color', color_order(2,:), 'LineWidth', 1, 'MarkerSize', 14); 
plot(flips(5:5:end)*180/pi, bssfp_approx_MT_noflow(2,5:5:end)  , '.', 'Color', color_order(3,:), 'LineWidth', 1, 'MarkerSize', 14); 
plot(flips(5:5:end)*180/pi, bssfp_approx_MT_flow(2,5:5:end)    , '.', 'Color', color_order(4,:), 'LineWidth', 1, 'MarkerSize', 14); 

xlabel('Flip angle (degrees)', 'FontSize', FontSize);
ylabel('Signal / M_0', 'FontSize', FontSize);
set(gca, 'FontSize', FontSize);
ylim([0 0.22]);
%legend('Exact  , no MT, no flow' , 'Exact  , no MT,      flow' , 'Exact  ,      MT, no flow' , 'Exact  ,      MT,      flow', ...
%       'Approx, no MT, no flow', 'Approx, no MT,      flow', 'Approx,      MT, no flow', 'Approx,      MT,      flow', 'location', 'southeast');
title('bSSFP', 'FontSize', FontSize);

subplot(1,2,2); hold on; set(gca, 'Box', 'On'); grid on;
plot(flips*180/pi, spgr_exact_noMT_noflow(2,:), '-', 'Color', color_order(1,:), 'LineWidth', 1); 
plot(flips*180/pi, spgr_exact_noMT_flow(2,:)  , '-', 'Color', color_order(2,:), 'LineWidth', 1); 
plot(flips*180/pi, spgr_exact_MT_noflow(2,:)  , '-', 'Color', color_order(3,:), 'LineWidth', 1); 
plot(flips*180/pi, spgr_exact_MT_flow(2,:)    , '-', 'Color', color_order(4,:), 'LineWidth', 1); 

plot(flips(5:5:end)*180/pi, spgr_approx_noMT_noflow(2,5:5:end), '.', 'Color', color_order(1,:), 'LineWidth', 1, 'MarkerSize', 14); 
plot(flips(5:5:end)*180/pi, spgr_approx_noMT_flow(2,5:5:end)  , '.', 'Color', color_order(2,:), 'LineWidth', 1, 'MarkerSize', 14); 
plot(flips(5:5:end)*180/pi, spgr_approx_MT_noflow(2,5:5:end)  , '.', 'Color', color_order(3,:), 'LineWidth', 1, 'MarkerSize', 14); 
plot(flips(5:5:end)*180/pi, spgr_approx_MT_flow(2,5:5:end)    , '.', 'Color', color_order(4,:), 'LineWidth', 1, 'MarkerSize', 14); 

xlabel('Flip angle (degrees)', 'FontSize', FontSize);
ylabel('Signal / M_0', 'FontSize', FontSize);
set(gca, 'FontSize', FontSize);
ylim([0 0.07]);
hLegend = legend('Exact  , no MT, no flow' , 'Exact  , no MT,      flow' , 'Exact  ,      MT, no flow' , 'Exact  ,      MT,      flow', ...
       'Approx, no MT, no flow', 'Approx, no MT,      flow', 'Approx,      MT, no flow', 'Approx,      MT,      flow', 'location', 'northeast');
set(hLegend, 'FontSize', 12);
title('SPGR', 'FontSize', FontSize);
%%
export_fig('steady_state_BMF_bSSFP_SPGR', '-r864', '-tif');
