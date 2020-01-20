function [mx,my,mz] = bloch_flow(b1, gr, taus, zetas, T1, T2, df, dp, mode, mx0, my0, mz0, F, lambda, m0, ASL_bolus)
% bloch_flow -- Bloch-Flow simulator using the extended matrix formalism
%  Usage
%    [mx,my,mz] = bloch_flow(b1, gr, taus, zetas, T1, T2, df, dp, mode, mx0, my0, mz0, F, lambda, m0, ASL_bolus)
%  Inputs
%    b1          RF pulse in G. Can be complex          [G]    (M x 1)
%    gr          1,2, or 3-dimensional gradient         [G/cm] (M x 1,2,or 3)
%    taus        time duration of each b1 and gr point  [sec]  (M x 1)
%    zetas       measurement time of each interval      [sec]  (M x 1)
%    T1          longitudinal relaxation time           [sec]  (1 x 1)
%    T2          transverse relaxation time             [sec]  (1 x 1)
%    df          array of off-resonance frequencies     [Hz]   (N x 1)
%    dp          array of spatial positions             [cm]   (P x 1,2,or 3)
%                Width should match width of gr.
%    mode        Bitmask mode:
%                Bit 0: 0-Simulate from start or M0, 1-Steady State (not supported)
%                Bit 1: 1-Record m at time points. 0-just end time.
%    mx0         array of starting x magnetization             (P x N)
%    my0         array of starting y magnetization             (P x N)
%    mz0         array of starting z magnetization             (P x N)
%    F           blood flow                             [mL/g/min]                   (1 x 1)
%    lambda      tissue/blood partition coefficient     [mL blood/g tissue]          (1 x 1)
%    m0          equilibrium magnetization              [magnetization/g tissue]     (1 x 1)
%    ASL_bolus   ASL bolus signal                       [magnetization/g tissue/sec] (M x 1)
%  Outputs
%    mx          array of the resulting x magnetization [magnetization/g tissue] (M x P x N) or (P x N)
%    my          array of the resulting y magnetization [magnetization/g tissue] (M x P x N) or (P x N)
%    mz          array of the resulting z magnetization [magnetization/g tissue] (M x P x N) or (P x N)
%
% Copyright (c) 2019. Namgyun Lee


%% Define function handles
I = eye(3); % 3 x 3 identity matrix

%--------------------------------------------------------------------------
% Relaxation operator
% Describes T1 and T2 relaxation over a time interval of duration "tau" [sec]
%--------------------------------------------------------------------------
E = @(tau) [exp(-tau/T2)        0              0     ;
                  0       exp(-tau/T2)         0     ;
                  0             0       exp(-tau/T1)];

%--------------------------------------------------------------------------
% Clearance operator
% Describes the clearance of transverse and longitudinal magnetization by
% venous flow over a time interval of duration "tau" [sec]
% F: [mL/g/min] * [min/60 sec] => [mL/g/sec]
%--------------------------------------------------------------------------
C = @(tau) [exp(-F/60*tau/lambda)            0                      0           ;
                     0             exp(-F/60*tau/lambda)            0           ;
                     0                       0            exp(-F/60*tau/lambda)];

%--------------------------------------------------------------------------
% Calculate the apparent T1
%--------------------------------------------------------------------------
T1app = 1 / (1 / T1 + (F / lambda / 60)); % apparent T1 [sec]

%% Perform Bloch-Flow simulation
gam = 4257.746778 * 2 * pi * 1e-2; % [rad/sec/uT]

M = length(b1); % number of time intervals
N = size(df,1); % number of off-resonance frequencies
P = size(dp,1); % number of spatial positions

m_zeta = zeros(M,3, 'double');

if mode == 2 % Record m at time points
    mx = zeros(M,P,N, 'double'); % M x P x N
    my = zeros(M,P,N, 'double'); % M x P x N
    mz = zeros(M,P,N, 'double'); % M x P x N
elseif mode == 0 % just end time
    mx = zeros(P,N, 'double'); % P x N
    my = zeros(P,N, 'double'); % P x N
    mz = zeros(P,N, 'double'); % P x N
end

for idx3 = 1:N
    df_ = df(idx3,1); % current off-resonance frequency [Hz]

    for idx2 = 1:P
        start_time = tic;
        fprintf('(%d/%d,%d/%d) Performing Bloch-Flow simulation ...', idx3, N, idx2, P);
        %------------------------------------------------------------------
        % Set up the starting magnetization
        %------------------------------------------------------------------
        ma = [mx0(idx2,idx3) my0(idx2,idx3) mz0(idx2,idx3)].';
        dp_ = dp(idx2,:);

        for idx1 = 1:M
            %==============================================================
            % Timing diagram
            %--------------------------------------------------------------
            %      b1(1)         b1(2)         b1(M-1)       b1(M)
            %        |             |             |             |
            %     ma | mb mc    md |             |             |
            %       \|/   |       \|             |             |
            %   -----+----x--------+----x--------+----x--------+----x-----
            %        t(1)          t(2)          t(M-1)        t(M)
            %        <------------><------------><------------><--------->
            %            tau(1)        tau(2)        tau(M-1)      tau(M)
            %        <--->         <--->         <--->         <--->
            %        zeta(1)       zeta(2)       zeta(M-1)     zeta(M)
            %==============================================================
            tau  = taus(idx1);  % current duration of the interval
            zeta = zetas(idx1); % current measurement time

            %--------------------------------------------------------------
            % RF excitation
            %--------------------------------------------------------------
            % Use Rodrigues' rotation formula
            b1x = real(b1(idx1)); % [G]
            b1y = imag(b1(idx1)); % [G]
            Bz = (gr(idx1,:) * dp_.') + 2 * pi * df_ / (gam * 1e2); % [G]
            B_abs = sqrt(b1x^2 + b1y^2 + Bz^2); % [G]
            theta = gam * 1e2 * B_abs * tau; % [rad/sec/uT] * [1e2uT/G] * [G] => [rad/sec]

            if B_abs > 0
                ux = b1x / B_abs;
                uy = b1y / B_abs;
                uz = Bz  / B_abs;
            else
                ux = 0;
                uy = 0;
                uz = 0;
            end

            % Rotations are left-handed
            c = cos(-theta);
            s = sin(-theta);
            one_minus_c = 1 - c;

            R11  = c + ux * ux * one_minus_c;
            R22  = c + uy * uy * one_minus_c;
            R33  = c + uz * uz * one_minus_c;
            R12a = ux * uy * one_minus_c;
            R12b = uz * s;
            R13a = ux * uz * one_minus_c;
            R13b = uy * s;
            R23a = uy * uz * one_minus_c;
            R23b = ux * s;

            R = [R11        , R12a - R12b, R13a + R13b ;
                 R12a + R12b, R22        , R23a - R23b ;
                 R13a - R13b, R23a + R23b, R33        ];

            %--------------------------------------------------------------
            % Compare Rodrigues' rotation formula with the matrix exponential
            %--------------------------------------------------------------
            %Omega = gam * 1e2 * [0 Bz -b1y; -Bz 0 b1x; b1y -b1x 0]; % [rad/sec/uT] * [1e2uT/G] * [G] => [rad/sec]
            %R2 = expm(Omega * tau); % [rad/sec] * [sec] => [rad]

            mb = R * ma;

            %--------------------------------------------------------------
            % Labeled flow
            % Calculate a new pseudo M0 term: M0 + s(t) * T1app, which is
            % the hypothetical equilibrium magnetization for modeling flow
            % during this short interval.
            %--------------------------------------------------------------
            pseudo_M0 = [0 0 m0 + ASL_bolus(idx1) * T1app].'; % 3 x 1

            %--------------------------------------------------------------
            % Relaxation (E) - Flow (F)
            %--------------------------------------------------------------
            CE = C(zeta) * E(zeta);
            mc = CE * mb + (I - CE) * pseudo_M0;

            %--------------------------------------------------------------
            % Record the magnetization at zeta
            %--------------------------------------------------------------
            m_zeta(idx1,:) = mc;

            %--------------------------------------------------------------
            % Relaxation (E) - Flow (F)
            %--------------------------------------------------------------
            CE = C(tau - zeta) * E(tau - zeta);
            md = CE * mc + (I - CE) * pseudo_M0;
            ma = md;
        end

        if mode == 2 % Record m at time points
            mx(:,idx2,idx3) = m_zeta(:,1); % M x P x N
            my(:,idx2,idx3) = m_zeta(:,2); % M x P x N
            mz(:,idx2,idx3) = m_zeta(:,3); % M x P x N
        elseif mode == 0 % just end time
            mx(idx2,idx3) = m_zeta(M,1); % P x N
            my(idx2,idx3) = m_zeta(M,2); % P x N
            mz(idx2,idx3) = m_zeta(M,3); % P x N
        end            
        fprintf('done! (%5.3f sec)\n' , toc(start_time));
    end
end

end