% demo_figure3_PCASL_low_parameter_sweep.m
% Written by Namgyun Lee
% Email: ggang56@gmail.com
% Started: 09/02/2018, Last modified: 12/21/2019

%% Clean slate
close all; clear all; clc;

%% Add paths
addpath(genpath('D:\ASL_project\mfiles_nam\Bloch_Flow_MT'));

%% Define parameters
% Used the recommended parameters in the consensus paper
F      = 0.3;         % brain blood flow [mL/g/min]
T1b    = 1650e-3;     % longitudinal relaxation time of arterial blood at 3T [sec]
T1     = 1200e-3;     % spin-lattice relaxation time of human brain at 3T [sec]
T2     = 1e10;        % spin-spin relaxation time of human brain at 3T [sec]
lambda = 0.9;         % brain/blood partition coefficient [mL of blood /g of brain tissue]
alpha  = 0.85;        % labeling efficiency for PCASL
alpha0 = 2 * alpha;   % labeling efficiency; 1 for saturation, 2 for inversion, and 0 for control
m0     = 1;           % equilibrium z magnetization
m0b    = m0 / lambda; % magnetization of blood
df     = 0;           % tissue off-resonance frequency [Hz]

T = 4;                % end time of a signal evolution [sec]
t_labeling = 0;       % start time of labeling RF pulses [sec]
nr_labeling = length(t_labeling);
labeling_on = true(nr_labeling,1);

%% Parameter sweep
F_range  = F; % [mL/g/min]
TD_range = (500:1:1500)*1e-3.';   % [sec]
TW_range = (1500:10:2000)*1e-3.'; % [sec]
T1_range = T1;
dt_range = (1:0.1:50)*1e-3.'; % sampling interval [sec]

F_length  = length(F_range);
TD_length = length(TD_range);
TW_length = length(TW_range);
T1_length = length(T1_range);
dt_length = length(dt_range);
nr_elements = F_length * TD_length * TW_length * T1_length;
space_size = [F_length, TD_length, TW_length, T1_length];

F_digits  = floor(log10(F_length))  + 1;
TD_digits = floor(log10(TD_length)) + 1;
TW_digits = floor(log10(TD_length)) + 1;
T1_digits = floor(log10(T1_length)) + 1;
dt_digits = floor(log10(dt_length)) + 1;

space            = zeros(nr_elements,5,dt_length, 'double');
nrmse            = zeros(nr_elements,dt_length, 'double');
max_deviation    = zeros(nr_elements,dt_length, 'double');
computation_time = zeros(nr_elements,dt_length, 'double');

for k = 1:dt_length
    start_time = tic;
    fprintf('Processing (%d/%d) ...', k, dt_length);
    dt = dt_range(k);    % sampling interval [sec]
    nt = floor(T / dt);  % number of samples
    t = (0:nt-1).' * dt;
    b1 = complex(zeros(nt,1, 'double'));
    gr = zeros(nt,1, 'double');
    taus  = dt * ones(nt,1, 'double');
    zetas = taus;

    for idx = 1:nr_elements
        %fprintf('Processing (%d/%d), (%d/%d) ...', k, dt_length, idx, nr_elements);
        [s1,s2,s3,s4] = ind2sub(space_size, idx);
        F  = F_range(s1);
        TD = TD_range(s2);
        TW = TW_range(s3);
        T1 = T1_range(s4);
        space(idx,:,k) = [F TD TW T1 dt];

        %------------------------------------------------------------------
        % Calculate the continuous ASL signal using the GKM
        %------------------------------------------------------------------
        PCASL_analytic = calculate_GKM_continuous_ASL_signal(t, F, lambda, m0b, alpha, T1, T1b, t_labeling, TD, TW);

        %------------------------------------------------------------------
        % Calculate the ASL bolus signal s(t)
        %------------------------------------------------------------------
        continuous_ASL_bolus = calculate_continuous_ASL_bolus(t, F, lambda, m0, alpha0, T1b, t_labeling, TD, TW * ones(nr_labeling,1), labeling_on); % [M0 / msec]

        %------------------------------------------------------------------
        % Perform Bloch simulation for the Bloch-Flow equations
        %------------------------------------------------------------------
        computation_start_time = tic;
        [mx1,my1,mz1] = bloch_flow_mex(b1, gr, taus, zetas, T1, T2, 0, 0, 2, 0, 0, 1, F, lambda, m0, zeros(nt,1));
        [mx2,my2,mz2] = bloch_flow_mex(b1, gr, taus, zetas, T1, T2, 0, 0, 2, 0, 0, 1, F, lambda, m0, continuous_ASL_bolus);
        PCASL_numeric = mz1 - mz2; % control - labeled

        %------------------------------------------------------------------
        % Calculate the computation time
        %------------------------------------------------------------------
        computation_time(idx,k) = toc(computation_start_time);

        %------------------------------------------------------------------
        % Calculate the maximum deviation
        %------------------------------------------------------------------
        max_deviation(idx,k) = max(PCASL_analytic - PCASL_numeric);

        %------------------------------------------------------------------
        % Calculate the normalized root-mean-square error
        %------------------------------------------------------------------
        nrmse(idx,k) = norm(PCASL_analytic - PCASL_numeric) / norm(PCASL_analytic);
        %fprintf('done! (%5.3f sec)\n', toc(start_time));
    end
    fprintf('done! (%5.3f sec)\n', toc(start_time));
end

%% Save results
save('PCASL_low_parameter_sweep', ...
    'F_range', 'TD_range', 'TW_range', 'T1_range', 'dt_range', 'nrmse', 'max_deviation', 'computation_time', '-v7.3');
