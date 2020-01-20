function [mxf,myf,mzf,mzs] = bloch_flow_MT(b1, gr, taus, zetas, T1f, T2f, df, dp, mode, mx0f, my0f, mz0f, mz0s, F, lambda, M0f, ASL_bolus, T1s, f, k, W)
% bloch_flow_MT -- Bloch-McConnell-Flow simulator using the extended matrix formalism
%  Usage
%    [mxf,myf,mzf,mzs] = bloch_flow_MT(b1, gr, taus, zetas, T1f, T2f, df, dp, mode, mx0f, my0f, mz0f, mz0s, F, lambda, M0f, ASL_bolus, T1s, f, k, W)
%  Inputs
%    b1          RF pulse in G. Can be complex          [G]       (M x 1)
%    gr          1,2, or 3-dimensional gradient         [G/cm]    (M x 1,2,or 3)
%    taus        time duration of each b1 and gr point  [sec]     (M x 1)
%    zetas       measurement time of each interval      [sec]     (M x 1)
%    T1f         free pool longitudinal relaxation time [sec]     (1 x 1)
%    T2f         free pool transverse relaxation time   [sec]     (1 x 1)
%    df          array of off-resonance frequencies     [Hz]      (N x 1)
%    dp          array of spatial positions             [cm]      (P x 1,2,or 3)
%                Width should match width of gr.
%    mode        Bitmask mode:
%                Bit 0: 0-Simulate from start or M0, 1-Steady State (not supported)
%                Bit 1: 1-Record m at time points. 0-just end time.
%    mx0f        array of starting x magnetization of the free pool      (P x N)
%    my0f        array of starting y magnetization of the free pool      (P x N)
%    mz0f        array of starting z magnetization of the free pool      (P x N)
%    mz0s        array of starting z magnetization of the semisolid pool (P x N)
%    F           blood flow                                         [mL blood/g tissue/min]      (1 x 1)
%    lambda      tissue/blood partition coefficient                 [mL blood/g tissue]          (1 x 1)
%    M0f         equilibrium magnetiztion of the free pool          [magnetization/g tissue]     (1 x 1)
%    ASL_bolus   ASL bolus signal                                   [magnetization/g tissue/sec] (M x 1)
%    T1s         semisolid longitudinal relaxation time             [sec]                        (1 x 1)
%    f           the fractional size of the semisolid pool                                       (1 x 1)
%    k           fundamental rate constant between the pools        [1/sec]                      (1 x 1)
%    G           array of the absorption lineshape of the semisolid [sec]                        (M x 1)
%    W           array of time-dependent saturation rates           [rad/sec]                    (M x 1)
%                W(t) = pi * (gam * 1e2)^2 * (b1x(t)^2 + b1y(t)^2) * G(t)
%                            ([rad/sec/uT] * [1e2uT/G])^2 * [G]^2 * [sec]
%  Outputs
%    mxf         array of the resulting x magnetization of the free pool      [magnetization/g tissue] (M x P x N) or (P x N)
%    myf         array of the resulting y magnetization of the free pool      [magnetization/g tissue] (M x P x N) or (P x N)
%    mzf         array of the resulting z magnetization of the free pool      [magnetization/g tissue] (M x P x N) or (P x N)
%    mzs         array of the resulting z magnetization of the semisolid pool [magnetization/g tissue] (M x P x N) or (P x N)
%
% Copyright (c) 2019. Namgyun Lee

%% Unpack tissue parameters
M0s = f * M0f;  % equilibrium magnetiztion of the semisolid pool
kf  = k * M0s;  % pseudo-first-order rate constant (free->semi) [1/sec]
ks  = k * M0f;  % pseudo-first-order rate constant (semi->free) [1/sec]

%% Define evolution operators
I = eye(4); % 4 x 4 identity matrix

%--------------------------------------------------------------------------
% Calculate the exchange operator
%--------------------------------------------------------------------------
A = @(tau) [f + 1     0                   0                            0                ;
              0     f + 1                 0                            0                ;
              0       0     1 + f * exp(-(f + 1) * ks * tau)   1 - exp(-(f + 1) * ks * tau) ;
              0       0     f - f * exp(-(f + 1) * ks * tau)   f + exp(-(f + 1) * ks * tau)] * 1 / (f + 1);
                        
                        
%--------------------------------------------------------------------------
% Relaxation operator
% Describes T1 and T2 relaxation over a time interval of duration "tau" [sec]
%--------------------------------------------------------------------------
E = @(tau) [exp(-tau / T2f)           0                  0                  0        ;
                   0           exp(-tau / T2f)           0                  0        ;
                   0                  0           exp(-tau / T1f)           0        ;
                   0                  0                  0           exp(-tau / T1s)];

%--------------------------------------------------------------------------
% Clearance operator
% Describes the clearance of transverse and longitudinal magnetization by
% venous flow over a time interval of duration "tau" [sec]
% F: [mL/g/min] * [min/60 sec] => [mL/g/sec]
%--------------------------------------------------------------------------
C = @(tau) [exp(-F/60*tau/lambda)              0                         0             0 ;
                      0              exp(-F/60*tau/lambda)               0             0 ;
                      0                        0              exp(-F/60*tau/lambda)    0 ;
                      0                        0                         0             1];

%--------------------------------------------------------------------------
% Calculate the apparent T1
% F / lambda => [mL blood/g tissue/min] / [mL blood/g tissue] * [min/60sec]
% => [1/sec]
%--------------------------------------------------------------------------
T1app = 1 / (1 / T1f + (F / lambda / 60)); % apparent T1 [sec]

%% Perform Bloch-Flow simulation
gam = 4257.746778 * 2 * pi * 1e-2; % [rad/sec/uT]

M = length(b1); % number of time intervals
N = size(df,1); % number of off-resonance frequencies
P = size(dp,1); % number of spatial positions

m_zeta = zeros(M,4, 'double');

if mode == 2 % Record m at time points
    mxf = zeros(M,P,N, 'double'); % M x P x N
    myf = zeros(M,P,N, 'double'); % M x P x N
    mzf = zeros(M,P,N, 'double'); % M x P x N
    mzs = zeros(M,P,N, 'double'); % M x P x N
elseif mode == 0 % just end time
    mxf = zeros(P,N, 'double'); % P x N
    myf = zeros(P,N, 'double'); % P x N
    mzf = zeros(P,N, 'double'); % P x N
    mzs = zeros(P,N, 'double'); % P x N
end

for idx3 = 1:N
    df_ = df(idx3,1); % current off-resonance frequency [Hz]

    for idx2 = 1:P
        start_time = tic;
        fprintf('(%d/%d,%d/%d) Performing Bloch-McConnell-Flow simulation (expm) ...', idx3, N, idx2, P);
        %------------------------------------------------------------------
        % Set up the starting magnetization (4 x 1)
        %------------------------------------------------------------------
        ma = [mx0f(idx2,idx3) my0f(idx2,idx3) mz0f(idx2,idx3) mz0s(idx2,idx3)].';
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
            % R: rotation for the free water and a saturation term for the
            % semisolid pool
            %--------------------------------------------------------------
            % Use Rodrigues' rotation formula
            b1x = real(b1(idx1)); % [G]
            b1y = imag(b1(idx1)); % [G]
            Bz = (gr(idx1,:) * dp_.') + 2 * pi * df_ / (gam * 1e2); % [G]
            B_abs = sqrt(b1x^2 + b1y^2 + Bz^2); % [G]

            % Calculate the rotation axis and angle
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

            R = [R11        , R12a - R12b, R13a + R13b, 0             ;
                 R12a + R12b, R22        , R23a - R23b, 0             ;
                 R13a - R13b, R23a + R23b, R33        , 0             ;
                 0          , 0          , 0          , exp(-W(idx1) * tau)];

            mb = R * ma;

            %--------------------------------------------------------------
            % Calculate a new pseudo M0 term:
            % -inv(Lambda + Gamma + Xi) * D(t)
            %
            % Lambda = [0 0   0   0 , Gamma = [-F/lambda     0        0     0
            %           0 0   0   0                0     -F/lambda    0     0
            %           0 0 -kf  ks                0         0     F/lambda 0
            %           0 0  kf -ks]               0         0        0     0]
            %
            % Xi = [-1/T2f    0      0      0     D(t) = [   0
            %          0   -1/T2f    0      0                0
            %          0      0   -1/T1f    0             (1/T1f + F/lambda)*M0f + s(t)
            %          0      0      0   -1/T1s]           1/T1s*M0s]
            %
            %--------------------------------------------------------------
            c = (1 + T1app * kf + T1s * ks); % denominator
            pseudo_M0 = [0;
                         0;
                         (1 + T1s * ks) / c * (M0f + ASL_bolus(idx1) * T1app) +      T1app * ks  / c * M0s; ...
                              T1s * kf  / c * (M0f + ASL_bolus(idx1) * T1app) + (1 + T1app * kf) / c * M0s];

            %--------------------------------------------------------------
            % Exchange (A) - Flow (C) - Relaxation (E)
            % expm((Lambda + Gamma + Xi) * tau)
            % ~= expm(Lambda * tau) * expm(Gamma * tau) * expm(Xi * tau)
            %  = A(tau) * C(tau) * E(tau)
            %--------------------------------------------------------------
            %ACE = A(zeta) * C(zeta) * E(zeta);
            Lambda = [0 0 0 0; 0 0 0 0; 0 0 -kf ks; 0 0 kf -ks]; % [1/sec]
            Gamma  = [-F/lambda/60 0 0 0; 0 -F/lambda/60 0 0; 0 0 -F/lambda/60 0; 0 0 0 0]; % [1/sec]
            Xi     = [-1/T2f 0 0 0; 0 -1/T2f 0 0; 0 0 -1/T1f 0; 0 0 0 -1/T1s]; % [1/sec]
            
            ACE = expm((Lambda + Gamma + Xi) * zeta);
            mc = ACE * mb + (I - ACE) * pseudo_M0;

            %--------------------------------------------------------------
            % Record the magnetization at zeta
            %--------------------------------------------------------------
            m_zeta(idx1,:) = mc;

            %--------------------------------------------------------------
            % Exchange (A) - Flow (C) - Relaxation (E)
            %--------------------------------------------------------------
            %ACE = A(tau - zeta) * C(tau - zeta) * E(tau - zeta);
            ACE = expm((Lambda + Gamma + Xi) * (tau - zeta));
            md = ACE * mc + (I - ACE) * pseudo_M0;
            ma = md;
        end

        if mode == 2 % Record m at time points
            mxf(:,idx2,idx3) = m_zeta(:,1); % M x P x N
            myf(:,idx2,idx3) = m_zeta(:,2); % M x P x N
            mzf(:,idx2,idx3) = m_zeta(:,3); % M x P x N
            mzs(:,idx2,idx3) = m_zeta(:,4); % M x P x N
        elseif mode == 0 % just end time
            mxf(idx2,idx3) = m_zeta(M,1); % P x N
            myf(idx2,idx3) = m_zeta(M,2); % P x N
            mzf(idx2,idx3) = m_zeta(M,3); % P x N
            mzs(idx2,idx3) = m_zeta(M,4); % P x N
        end
        fprintf('done! (%5.3f sec)\n' , toc(start_time));
    end
end

end