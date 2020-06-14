% demo_Su_2017_MRM_MRF_ASL.m
% Written by Namgyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 10/07/2019, Last modified: 06/13/2020

%% Clean slate
close all; clear all; clc;

%% Add paths
addpath(genpath('..\..\Bloch_Flow_MT'));

%% Initialize a random number generator
rng('default');

%% Define simulation parameters
F      = 0.6;         % cerebral blood flow [mL/g/min] WM: 34 [mL/100 mL/min], GM: 42 [mL/100 mL/min]
T1b    = 1664e-3;     % longitudinal relaxation time of arterial blood at 3T [sec]
T1     = 1200e-3;     % spin-lattice relaxation time of human brain at 3T [sec] WM: 1110 msec, GM: 1470 msec
T2     = 80e-3;       % spin-spin relaxation time of human brain at 3T [sec]
lambda = 0.9;         % brain/blood partition coefficient [mL of blood /g of brain tissue]
TD     = 1000e-3;     % transit delay or tissue bolus arrival time (BAT) (deltat) [sec]
M0     = 1;           % equilibrium z magnetization [magnetization/g tissue]
M0b    = M0 / lambda; % M0 arterial blood = M0 tissue / lambda
df     = 30;          % tissue off-resonance frequency [Hz]

% Assume pulsed ASL for MRF-ASL
alpha  = 1;           % labeling efficiency for PASL
alpha0 = 2 * alpha;   % labelling efficiency; 1 for saturation, 2 for inversion, and 0 for control
B1     = 1;           % transmit sensitivity
flip   = 40 * pi / 180; % nominal flip angle [rad]
dt     = 1e-3;        % sampling interval in time [sec]

%% Calculate the labeling duration
nr_TRs = 30; % number of TRs
TWs = floor(189 * cos(2 * pi * (0:nr_TRs-1).' / (2 * nr_TRs)) + 261) * 1e-3; % labeling duration [sec]

%% Calculate TDs [sec]
TDs = TD * ones(nr_TRs,1, 'double'); % [sec]

%% Calculate TRs [sec]
TRs = TWs;

%% Calculate the label index
nr_labels = nr_TRs / 2;
random_index = randperm(nr_TRs).';
label_index = sort(random_index(1:nr_labels));

%% Calculate the start time of a labeling (Label and Control) RF pulse [sec]
t_labeling = zeros(nr_TRs,1, 'double'); % [sec]
for idx = 1:nr_TRs
    t_labeling(idx) = sum(TRs(1:idx-1)); 
end

labeling_on = false(nr_TRs,1);
labeling_on(label_index) = true;

%% Calculate sampling times [sec]
N = sum(ceil(TRs / dt)); % the number of time intervals
t = zeros(N,1, 'double');

index_offset = 0;
for idx = 1:nr_TRs
    nr_samples = ceil(TRs(idx) / dt);
    index_range = (1:nr_samples).' + index_offset;
    t(index_range) = (0:nr_samples-1).' * dt + t_labeling(idx);
    index_offset = index_offset + nr_samples;
end

%% Calculate the ASL signal
ASL_signal = zeros(N,1, 'double');

%--------------------------------------------------------------------------
% Calculate the apparent T1
%--------------------------------------------------------------------------
T1app = 1 / (1 / T1 + (F / lambda / 60)); % apparent T1 [sec]

%--------------------------------------------------------------------------
% Calculate k [1/sec]
%--------------------------------------------------------------------------
k = 1 / T1b - 1 / T1app; % [1/sec]

%--------------------------------------------------------------------------
% Add all Label pulses to the current time
%--------------------------------------------------------------------------
for idx2 = 1:N
    t1 = t(idx2); % current time [sec]
    for idx1 = 1:nr_TRs
        if labeling_on(idx1)
            t_label = t_labeling(idx1);
            TD = TDs(idx1);
            TW = TWs(idx1);

            if (t1 <= t_label + TD)
                ASL_signal(idx2) = ASL_signal(idx2) + 0;

            elseif ((t1 > t_label + TD) && (t1 < t_label + TD + TW))
                % f [mL/g/min]: [mL/g/min] * [min/60 sec] => [mL/g/sec]
                ASL_signal(idx2) = ASL_signal(idx2) - (F / lambda / 60) * M0 * alpha0 * exp(-(t1 - t_label) / T1b) * ...
                    exp(k * (t1 - t_label)) * (exp(-k * TD) - exp(-k * (t1 - t_label))) / k;

            elseif (t1 >= t_label + TD + TW)
                ASL_signal(idx2) = ASL_signal(idx2) - (F / lambda / 60) * M0 * alpha0 * exp(-(t1 - t_label) / T1b) * ...
                    exp(k * (t1 - t_label)) * (exp(-k * TD) - exp(-k * (TD + TW))) / k;
            end
        end
    end
end

%% Solve the Bloch equations in a piecewise manner
Mtissue  = zeros(N,1, 'double');
Mtissue0 = zeros(nr_TRs,1, 'double');
Mxy      = zeros(nr_TRs,1, 'double'); % measured signal in MRF-ASL

Mtissue0_block = M0; % M0
index_offset = 0;

for idx = 1:nr_TRs
    %----------------------------------------------------------------------
    % Set the TR of a current block
    %----------------------------------------------------------------------
    TR = TRs(idx);

    %----------------------------------------------------------------------
    % Calculate the number of samples in a current block
    %----------------------------------------------------------------------
    nr_samples = ceil(TR / dt);

    %----------------------------------------------------------------------
    % current block start time [sec]
    %----------------------------------------------------------------------
    t_block_start = t_labeling(idx);

    %----------------------------------------------------------------------
    % current Mz index range
    %----------------------------------------------------------------------
    index_range = (1:nr_samples).' + index_offset;

    %----------------------------------------------------------------------
    % Calculate time samples within a current block
    %----------------------------------------------------------------------
    t_block = t(index_range);

    %----------------------------------------------------------------------
    % Record Mtissue at t=0 in a TR
    %----------------------------------------------------------------------
    Mtissue0(idx) = Mtissue0_block;

    %----------------------------------------------------------------------
    % Calculate the magnetization evolution within the current block
    %----------------------------------------------------------------------
    Mtissue_block = Mtissue0_block * exp(-(t_block - t_block_start) / T1app) + ...
                    M0 * (1 - exp(-(t_block - t_block_start) / T1app)) - ...
                    ASL_signal(index_range(1)) * exp(-(t_block - t_block_start) / T1app) + ASL_signal(index_range);

    %----------------------------------------------------------------------
    % Calculate the Mtissue at the end of current TR
    %----------------------------------------------------------------------
    if idx < nr_TRs
        t1 = t_block_start + TR;
        Mtissue_end = Mtissue0_block * exp(-(t1 - t_block_start) / T1app) + ...
                      M0 * (1 - exp(-(t1 - t_block_start) / T1app)) - ...
                      ASL_signal(index_range(1)) * exp(-(t1 - t_block_start) / T1app) + ASL_signal(index_range(end)+1);
    end

    %----------------------------------------------------------------------
    % Mtissue at t=0 in a TR is determined by Mtissue * cos(alpha) at the
    % end of the previous TR
    %----------------------------------------------------------------------
    Mtissue0_block = Mtissue_end * cos(flip * B1);

    %----------------------------------------------------------------------
    % The measured signal in MRF-ASL is M * sin(alpha)
    %----------------------------------------------------------------------
    Mxy(idx) = Mtissue_end * sin(flip * B1);

    %----------------------------------------------------------------------
    % Record the magnetization of the tissue compartment
    %----------------------------------------------------------------------
    Mtissue(index_range) = Mtissue_block;
    index_offset = index_offset + nr_samples;
end

%% Calculate the ASL bolus signal
ASL_bolus = calculate_pulsed_ASL_bolus(t, F, lambda, M0, alpha0, T1b, t_labeling, TDs, TWs, labeling_on); % [magnetization/g tissue/sec]

%% Calculate b1 [uT], gradient [G/cm], and taus [sec]
gam = 4257.746778 * 2 * pi * 1e-2; % [rad/sec/uT]

b1   = zeros(N,1, 'double'); % [uT]
gr   = zeros(N,1, 'double'); % [G/cm]
taus = zeros(N,1, 'double'); % [sec]

index_offset = 0;
for idx = 1:nr_TRs
    %----------------------------------------------------------------------
    % Timing diagram:
    %
    % labeling1          labeling2
    % |                  |                  |
    % |     BLOCK1      RF1     BLOCK2     RF2
    % |                  |                  |
    % |        TR1       |        TR2       |
    % |<---------------->|<---------------->|
    % |                  |                  |
    % +------+------+------+------+------+------+------+------+------> time
    % 0      dt    2dt    3dt    4dt    5dt    6dt 
    % x      x      x    x      x      x      : 6 start times (t)
    % |<---->|<---->|<-->|<---->|<---->|<-->|
    %    1      2     3                       : 3 timesteps (tau)
    %                       1      2     3    : 3 timesteps (tau)
    % N: ceil(5.5dt/dt) = 6 timesteps
    %----------------------------------------------------------------------
    nr_timesteps = ceil(TRs(idx) / dt); % equal to the number of samples per a timestep
    index_range = (1:nr_timesteps).' + index_offset;

    %----------------------------------------------------------------------
    % Calculate b1
    %----------------------------------------------------------------------
    if idx ~= nr_TRs
        % [rad] / [sec] / [rad/sec/uT] => [uT]
        b1(index_range(end)+1) = (flip * B1) / dt / gam; % [uT]
    end

    %----------------------------------------------------------------------
    % Calculate the duration (tau)
    %----------------------------------------------------------------------
    taus(index_range) = cat(1, dt * ones(nr_timesteps-1,1, 'double'), TRs(idx) - (nr_timesteps-1) * dt); % [sec]

    %----------------------------------------------------------------------
    % Update the starting index of a next block
    %----------------------------------------------------------------------
    index_offset = index_offset + nr_timesteps;
end

%% Perform Bloch-Flow simulation
% Note: If T2 is incorrectly set up, then it is not possible to get the
% correct signal evolution!!
tic; [mx_numeric_without_T2,my_numeric_without_T2,mz_numeric_without_T2] = bloch_flow(b1*1e-2, gr, taus, 0*taus, T1, 1e-20, 0, 0, 2, 0, 0, 1, F, lambda, M0, ASL_bolus); toc;
mxy_numeric_without_T2 = mx_numeric_without_T2 + 1j*my_numeric_without_T2;

tic; [mx_numeric_with_T2,my_numeric_with_T2,mz_numeric_with_T2] = bloch_flow(b1*1e-2, gr, taus, 0*taus, T1, T2, 0, 0, 2, 0, 0, 1, F, lambda, M0, ASL_bolus); toc;
mxy_numeric_with_T2 = mx_numeric_with_T2 + 1j*my_numeric_with_T2;

tic; [mx_numeric_with_T2_df,my_numeric_with_T2_df,mz_numeric_with_T2_df] = bloch_flow(b1*1e-2, gr, taus, 0*taus, T1, T2, df, 0, 2, 0, 0, 1, F, lambda, M0, ASL_bolus); toc;
mxy_numeric_with_T2_df = mx_numeric_with_T2_df + 1j*my_numeric_with_T2_df;

%% Display label and control blocks
figure('Color', 'w', 'Position', [6 2 776 814]);
color_order1 = cbrewer('qual', 'Set1', 8);
color_order2 = get(gca, 'colororder');

label_color   = [120 171 48] / 255; % green
control_color = [255 255 255] / 255; % white

%--------------------------------------------------------------------------
% Display Labeling duration
%--------------------------------------------------------------------------
h1 = subplot(4,1,1); hold on; set(gca, 'Box', 'On'); grid on;
set(gca, 'Position', [0.13 0.805 0.775 0.176]);
index = (1:nr_TRs).';
plot(TWs*1e3, '.-', 'Color', 'k', 'LineWidth', 1, 'MarkerSize', 16, 'HandleVisibility', 'Off');
plot(index(labeling_on) , TWs(labeling_on)*1e3 , '.', 'Color', color_order1(1,:), 'MarkerSize', 16);
plot(index(~labeling_on), TWs(~labeling_on)*1e3, '.', 'Color', color_order2(1,:), 'MarkerSize', 16);
xlim([1 nr_TRs]);
ylim([0 500]);
legend('Label', 'Control');
xlabel('Image index');
ylabel('Duration (msec)');
title('Labeling Duration (Label + Control)');

%--------------------------------------------------------------------------
% Display label and control blocks
%--------------------------------------------------------------------------
h2 = subplot(4,1,2); hold on; set(gca, 'Box', 'On');
set(gca, 'Position', [0.13 0.55 0.775 0.176]); % [left bottom width height]
label_acquired   = true;
control_acquired = true;

for idx = 1:nr_TRs

    t_start = t_labeling(idx);
    t_end   = t_start + TRs(idx);

    x = [t_start t_start t_end t_end];
    y = [0 1.2 1.2 0];

    label_block = labeling_on(idx);

    if label_block
        face_color = label_color;
    else
        face_color = control_color;
    end

    h = patch('XData', x, 'YData', y, 'FaceColor', face_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'LineWidth', 1, 'LineStyle', '-');

    if label_block && label_acquired
        set(h, 'HandleVisibility', 'Off');
    end
    if ~label_block && control_acquired
        set(h, 'HandleVisibility', 'Off');
    end

    if label_block
        label_acquired = true;
    end
    if ~label_block
        control_acquired = true;
    end
end

%--------------------------------------------------------------------------
% Display Mtissue (GKM) vs Numeric 
%--------------------------------------------------------------------------
plot(t, mz_numeric_with_T2   , '-' , 'Color', 'k'              , 'LineWidth', 1);
plot(t, mz_numeric_with_T2_df, '-' , 'Color', color_order2(3,:), 'LineWidth', 1);
plot(t, Mtissue              , '-' , 'Color', color_order2(1,:), 'LineWidth', 1);
plot(t, mz_numeric_without_T2, '--', 'Color', color_order2(2,:), 'LineWidth', 1);
xlim([0 t(end)]);
ylim([0 1.2]);
xlabel('Time (sec)');
ylabel('Signal / M_0');
title('MRF Signal Evolution M_{z}^{tissue}');
legend_text = {sprintf('Numeric: M_z^{tissue} (T_2=%2.0f msec,\\Deltaf=%2.0f Hz)', T2*1e3, 0), sprintf('Numeric: M_z^{tissue} (T_2=%2.0f msec,\\Deltaf=%2.0f Hz)', T2*1e3, df), 'GKM: M_z^{tissue}', 'Numeric: M_z^{tissue} (T_2=0 msec,\Deltaf=0 Hz)'};
hLegend = legend(legend_text, 'Location', 'best', 'Orientation', 'vertical');
hLegend.NumColumns = 2;
set(hLegend, 'Position', [0.265 0.662 0.628 0.057]);

%--------------------------------------------------------------------------
% Display label and control blocks
%--------------------------------------------------------------------------
h3 = subplot(4,1,3); hold on; set(gca, 'Box', 'On');
set(gca, 'Position', [0.13 0.30 0.775 0.176]); % [left bottom width height]
label_acquired   = true;
control_acquired = true;

for idx = 1:nr_TRs

    t_start = t_labeling(idx);
    t_end   = t_start + TRs(idx);

    x = [t_start t_start t_end t_end];
    y = [0 1.2 1.2 0];

    label_block = labeling_on(idx);

    if label_block
        face_color = label_color;
    else
        face_color = control_color;
    end

    h = patch('XData', x, 'YData', y, 'FaceColor', face_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'LineWidth', 1, 'LineStyle', '-');

    if label_block && label_acquired
        set(h, 'HandleVisibility', 'Off');
    end
    if ~label_block && control_acquired
        set(h, 'HandleVisibility', 'Off');
    end

    if label_block
        label_acquired = true;
    end
    if ~label_block
        control_acquired = true;
    end
end

%--------------------------------------------------------------------------
% Display Mtissue (GKM) vs Numeric 
%--------------------------------------------------------------------------
plot(t, abs(mxy_numeric_with_T2)   , '-' , 'Color', 'k'              , 'LineWidth', 1);
plot(t, abs(mxy_numeric_with_T2_df), '-' , 'Color', color_order2(3,:), 'LineWidth', 1);
plot(t_labeling + TRs, Mxy         , '*' , 'Color', color_order2(1,:), 'LineWidth', 1);
plot(t, abs(mxy_numeric_without_T2), '--', 'Color', color_order2(2,:), 'LineWidth', 1);
xlim([0 t(end)]);
ylim([0 1.2]);
xlabel('Time (sec)');
ylabel('Signal / M_0');
title('MRF Signal Evolution M_{xy}^{tissue}');
legend_text = {sprintf('Numeric: M_{xy}^{tissue} (T_2=%2.0f msec,\\Deltaf=%2.0f Hz)', T2*1e3, 0), sprintf('Numeric: M_{xy}^{tissue} (T_2=%2.0f msec,\\Deltaf=%2.0f Hz)', T2*1e3, df), 'GKM: M_{xy}^{tissue}', 'Numeric: M_z^{tissue} (T_2=0 msec,\Deltaf=0 Hz)'};
hLegend = legend(legend_text, 'Location', 'northeast', 'Orientation', 'vertical');
hLegend.NumColumns = 2;

%--------------------------------------------------------------------------
% Display signal difference between Mtissue (GKM) vs Numeric 
%--------------------------------------------------------------------------
h4 = subplot(4,1,4); hold on; set(gca, 'Box', 'On');
set(gca, 'Position', [0.13 0.05 0.775 0.176]); % [left bottom width height]
plot(t,(Mtissue-mz_numeric_without_T2)*1e2, '-', 'LineWidth', 1); grid on;
xlim([0 t(end)]);
xlabel('Time (sec)');
ylabel('Signal / M_0 (%)');
title('Difference (GKM - Numeric (T_2=0 msec,\Deltaf=0 Hz))');

%% Save as a .tiff file
export_fig('../../figs/figure5', '-r600', '-tif');
