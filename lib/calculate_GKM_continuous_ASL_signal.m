function DELTAM = calculate_GKM_continuous_ASL_signal(t, f, lambda, M0b, alpha, T1, T1b, t_labeling, deltat, tau)
% calculate_GKM_continuous_ASL_signal -- Calculate the continuous labeling ASL signal
%  Usage
%    DELTAM = calculate_GKM_continuous_ASL_signal(t, f, lambda, M0b, alpha, T1, T1b, t_labeling, deltat, tau)
%  Inputs
%    t              vector of time points                          [sec] (N x 1)
%    f              blood flow                                     [mL blood/g tissue/min] => [mL/g/min]
%    lambda         tissue/blood partition coefficient for water   [mL blood/g tissue]
%    M0b            equilibrium magnetization of arterial blood
%    alpha          accounts for incomplete inversion during
%                   the tagging pulse
%    T1             longitudinal relaxation time of tissue water
%                   in the absence of flow                         [sec]
%    T1b            longitudinal relaxation time of arterial blood [sec]
%    t_labeling     vector of onset times of labeling              [sec] (M x 1)
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

% Calculate the continuous labeling ASL signal 
% Use Equation 5 in Buxton et al. 1998 MRM.
DELTAM = zeros(N,1, 'double');
q_ss   = zeros(N,1, 'double'); % steady-state

for idx2 = 1:nr_labeling
    ts = t_labeling(idx2);

    for idx1 = 1:N
        t1 = t(idx1); % current time [sec]

        if (t1 < ts + deltat)
            q_ss(idx1)    = 0;
            DELTAM(idx1) = DELTAM(idx1) + 0;

        elseif ((t1 >= ts + deltat) && (t1 < ts + tau + deltat))
            % f [mL/g/min]: [mL/g/min] * [min/60 sec] => [mL/g/sec]
            q_ss(idx1) = 1 - exp(-(t1 - ts - deltat) / T1app);
            DELTAM(idx1) = DELTAM(idx1) + 2 * M0b * (f / 60) * T1app * alpha * exp(-deltat / T1b) * q_ss(idx1);

        elseif (t1 >= ts + tau + deltat)
            q_ss(idx1) = 1 - exp(-tau / T1app);
            DELTAM(idx1) = DELTAM(idx1) + 2 * M0b * (f / 60) * T1app * alpha * exp(-deltat/ T1b) * exp(-(t1 - ts - tau - deltat) / T1app) * q_ss(idx1);
        end
    end
end

end