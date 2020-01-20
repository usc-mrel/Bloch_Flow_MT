function bolus = calculate_pulsed_ASL_bolus(t, F, lambda, m0, alpha0, T1b, t_labeling, TDs, TWs, labeling_on)
% calculate_pulsed_ASL_bolus -- Calculate the pulsed ASL bolus signal
%  Usage
%    bolus = calculate_pulsed_ASL_bolus(t, F, lambda, m0, alpha, T1b, t_labeling, TDs, TWs, labeling_on)
%  Inputs
%    t               vector of time points                          [sec]  (N x 1)
%    F               blood flow                                     [mL blood/g tissue/min] (e.g., 1 [mL/g/min])
%    lambda          tissue/blood partition coefficient for water   [mL blood/g tissue]
%    m0              equilibrium magnetization per g of tissue      [magnetization/g tissue]
%    alpha0          labelling efficiency; 1 for saturation, 2 for inversion, and 0 for control
%    T1b             longitudinal relaxation time of arterial blood [sec]
%    t_labeling      vector of onset times of RF labeling           [sec]  (M x 1)
%    TDs             vector of arterial transit times               [sec]  (M x 1)
%    TWs             vector of bolus widths of labeled blood        [sec]  (M x 1)
%    labeling_on     logical index for labeling; 0:no labeling, 1:labeling (M x 1)
%  Outputs
%    bolus           pulsed ASL bolus signal                        [magnetization/g tissue/sec] (N x 1)
%
% Copyright (c) 2019. Namgyun Lee

N = length(t);
nr_labeling = sum(labeling_on > 0);
labeling_index = find(labeling_on);

% Calculate the pulsed ASL bolus signal
bolus = zeros(N,1, 'double');
for idx = 1:nr_labeling

    % Get current parameters
    ts = t_labeling(labeling_index(idx)); % [sec]
    TD = TDs(labeling_index(idx));        % [sec]
    TW = TWs(labeling_index(idx));        % [sec]

    % F/lambda: [mL blood/g tissue/min] / [mL blood/g tissue] => [1/min]
    % [1/min] * [min/60 sec] => [1/sec]
    logical_index = (t >= ts + TD) & (t <= ts + TD + TW);
    bolus(logical_index) = bolus(logical_index) - (F / lambda / 60) * m0 * alpha0 * exp(-(t(logical_index) - ts) / T1b);
end

end