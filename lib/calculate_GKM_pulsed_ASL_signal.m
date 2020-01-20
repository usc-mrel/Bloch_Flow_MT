function DELTAM = calculate_GKM_pulsed_ASL_signal(t, f, lambda, M0b, alpha, T1, T1b, t_labeling, deltat, tau)
% calculate_GKM_pulsed_ASL_signal -- Calculate the pulsed labeling ASL signal
%  Usage
%    DELTAM = calculate_GKM_pulsed_ASL_signal(t, f, lambda, M0b, alpha, T1, T1b, t_labeling, deltat, tau)
%  Inputs
%    t              vector of time points                          [sec] (N x 1)
%    f              blood flow                                     [mL blood/g tissue/min] => [mL/g/min]
%    lambda         tissue/blood partition coefficient for water   [mL blood/g tissue]
%    M0b            equilibrium magnetization of arterial blood    [magnetization/mL blood]
%    alpha          accounts for incomplete inversion during
%                   the tagging pulse
%    T1             longitudinal relaxation time of tissue water
%                   in the absence of flow                         [sec]
%    T1b            longitudinal relaxation time of arterial blood [sec]
%    t_labeling     vector of onset times of RF labeling           [sec] (M x 1)
%    deltat         transit delay                                  [sec]
%    tau            bolus width (duration) of labeled blood        [sec]
%  Outputs
%    bolus          continuous ASL bolus signal                    [magnetization/g tissue/sec] (N x 1)
%
% Copyright (c) 2019. Namgyun Lee

N = length(t);
nr_labeling = length(t_labeling);

% Calculate the apparent relaxation time
T1app = 1 / (1 / T1 + (f / lambda / 60)); % apparent T1 [sec]

% Calculate k
k = 1 / T1b - 1 / T1app; % [1/sec]

% Calculate the pulsed labeling ASL signal 
% Use Equation 3 in Buxton et al. 1998 MRM.
DELTAM = zeros(N,1, 'double');
q_p    = zeros(N,1, 'double');

for idx2 = 1:nr_labeling
    ts = t_labeling(idx2);

    for idx1 = 1:N
        t1 = t(idx1); % current time [msec]

        if (t1 <= ts + deltat)
            q_p(idx1)    = 0;
            DELTAM(idx1) = DELTAM(idx1) + 0;

        elseif ((t1 > ts + deltat) && (t1 < ts + tau + deltat))
            % f [mL/g/min]: [mL/g/min] * [min/60 sec] => [mL/g/sec]
            q_p(idx1) = exp(k * (t1 - ts)) * (exp(-k * deltat) - exp(-k * (t1 - ts))) / (k * (t1 - ts - deltat));
            DELTAM(idx1) = DELTAM(idx1) + 2 * M0b * (f / 60) * (t1 - ts - deltat) * alpha * exp(-(t1 - ts) / T1b) * q_p(idx1);

        elseif (t1 >= ts + tau + deltat)
            q_p(idx1) = exp(k * (t1 - ts)) * (exp(-k * deltat) - exp(-k * (tau + deltat))) / (k * tau);
            DELTAM(idx1) = DELTAM(idx1) + 2 * M0b * (f / 60) * tau * alpha * exp(-(t1 - ts)/ T1b) * q_p(idx1);
        end
    end
end

end