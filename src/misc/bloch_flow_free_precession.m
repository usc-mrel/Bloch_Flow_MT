function [mx,my,mz] = bloch_flow_free_precession(b1, gr, taus, zetas, T1, T2, df, dp, mx0, my0, mz0, F, lambda, m0, ASL_bolus)
% bloch_flow -- Bloch-Flow simulator using the extended matrix formalism
%  Usage
%    [mx,my,mz] = bloch_flow(b1, gr, taus, zetas, T1, T2, df, dp, mx0, my0, mz0, F, lambda, m0, ASL_bolus)
%  Inputs
%    b1          RF pulse in G. Can be complex          [G]        (M x 1)
%    gr          1,2, or 3-dimensional gradient         [G/cm]     (M x 1,2,or 3)
%    taus        time duration of each b1 and gr point  [sec]      (M x 1)
%    zetas       measurement time of each interval      [sec]      (M x 1)
%    T1          longitudinal relaxation time           [sec]      (1 x 1)
%    T2          transverse relaxation time             [sec]      (1 x 1)
%    df          array of off-resonance frequencies     [Hz]       (N x 1)
%    dp          array of spatial positions             [cm]       (P x 1,2,or 3)
%                Width should match width of gr.
%    mx0         array of starting x magnetization                 (P x 1)
%    my0         array of starting y magnetization                 (P x 1)
%    mz0         array of starting z magnetization                 (P x 1)
%    F           blood flow                             [mL/g/min] (1 x 1)
%    lambda      tissue/blood partition coefficient     [mL blood/g tissue]          (1 x 1)
%    m0          equilibrium magnetization              [magnetization/g tissue]     (1 x 1)
%    ASL_bolus   ASL bolus signal                       [magnetization/g tissue/sec] (M x 1)
%  Outputs
%    mx          array of the resulting x magnetization [magnetization/g tissue] (M x P x N)
%    my          array of the resulting y magnetization [magnetization/g tissue] (M x P x N)
%    mz          array of the resulting z magnetization [magnetization/g tissue] (M x P x N)
%
% Copyright (c) 2019. Namgyun Lee


%% Define function handles
I = eye(3); % 3 x 3 identity matrix

%--------------------------------------------------------------------------
% Rotation matrices (Rotations are left-handed)
%--------------------------------------------------------------------------
% Rotation by alpha about the x-axis
Rx = @(alpha) [     1            0            0      ;
                    0       cos(alpha)    sin(alpha) ;
                    0      -sin(alpha)    cos(alpha)];

% Rotation by alpha about the y-axis
Ry = @(alpha) [ cos(alpha)      0        -sin(alpha) ;
                    0           1              0     ;
                sin(alpha)      0         cos(alpha)];

% Rotation by alpha about the z-axis
Rz = @(alpha) [ cos(alpha)   sin(alpha)        0     ;
               -sin(alpha)   cos(alpha)        0     ;
                    0           0              1    ];

%--------------------------------------------------------------------------                   
% Relaxation operator
% Describes T1 and T2 relaxation over a time interval of duration "tau" [sec]
%--------------------------------------------------------------------------                   
E = @(tau) [exp(-tau/T2)        0              0     ;
                  0       exp(-tau/T2)         0     ;
                  0             0       exp(-tau/T1)];

%--------------------------------------------------------------------------                   
% Free precession over a time interval of duration "tau" about the z-axis
% [2*pi rad/cycle] * [Hz] * [sec] => [rad]
%--------------------------------------------------------------------------                   
P = @(tau) [ cos(2*pi*df*tau)  sin(2*pi*df*tau)         0 ;
            -sin(2*pi*df*tau)  cos(2*pi*df*tau)         0 ;
                     0                0                 1];

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
% T1 Recovery + inflowing arterial flow of un-labeled blood
%--------------------------------------------------------------------------                   
D = @(tau) (I - C(tau) * E(tau)) * [0 0 m0].';

%% Perform Bloch-Flow simulation
gam = 4257.746778 * 2 * pi * 1e-2; % [rad/sec/uT]

M = length(b1); % number of time intervals
N = size(df,1); % number of off-resonance frequencies
P = size(dp,1); % number of spatial positions

m = zeros(3,M,P,N, 'double');

% b1 in G: [rad/sec/uT] * [1e2uT/G] * [G] * [sec] => [rad]
alpha = gam * 1e2 * abs(b1) .* taus; % RF pulse angle [rad]
phi   = angle(b1);                   % RF pulse phase [rad]

for idx3 = 1:N
    df_ = df(idx3,1); % current off-resonance frequency [Hz]

    for idx2 = 1:P
        start_time = tic;
        fprintf('(%d/%d,%d/%d) Performing Bloch-Flow simulation ...', idx3, N, idx2, P);
        %------------------------------------------------------------------
        % Set up the starting magnetization
        %------------------------------------------------------------------
        ma = [mx0(idx2,1) my0(idx2,1) mz0(idx2,1)].';
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
            %        <------------><------------><------------><------------>
            %            tau(1)        tau(2)        tau(M-1)      tau(M)
            %        <--->         <--->         <--->         <--->
            %        zeta(1)       zeta(2)       zeta(M-1)     zeta(M)
            %==============================================================
            tau  = taus(idx1);  % current duration of the interval
            zeta = zetas(idx1); % current measurement time

            %--------------------------------------------------------------
            % RF excitation
            %--------------------------------------------------------------
            mb = Rz(-phi(idx1)) * Rx(alpha(idx1)) * Rz(phi(idx1)) * ma;

            %--------------------------------------------------------------
            % Labeled flow
            % Calculate the amount of labeled inflow during the kth 
            % interval using the hard pulse approximation
            %--------------------------------------------------------------
            mb = mb + [0 0 ASL_bolus(idx1) * tau].'; % inflow of labeled blood

            %--------------------------------------------------------------
            % Relaxation (E) - Free precession (P) - Flow (F)
            %--------------------------------------------------------------
            % 1T = 1e4G = 1e6uT => G = 1e2uT
            % [rad/sec/uT] * [1e2uT/G] * [G/cm] * [cm] => [rad/sec]
            theta = (gam * 1e2 * (gr(idx1,:) * dp_.') + 2 * pi * df_) * zeta; % [rad]
            mc = C(zeta) * Rz(theta) * E(zeta) * mb + D(zeta);

            %--------------------------------------------------------------
            % Relaxation (E) - Free precession (P) - Flow (F)
            %--------------------------------------------------------------
            theta = (gam * 1e2 * (gr(idx1,:) * dp_.') + 2 * pi * df_) * (tau - zeta); % [rad]
            md = C(tau - zeta) * Rz(theta) * E(tau - zeta) * mc + D(tau - zeta);
            m(:,idx1,idx2,idx3) = mc;
            ma = md;
        end
        fprintf('done! (%5.3f sec)\n' , toc(start_time));
    end
end

mx = reshape(m(1,:,:,:), [M P N]); % M x P x N
my = reshape(m(2,:,:,:), [M P N]); % M x P x N
mz = reshape(m(3,:,:,:), [M P N]); % M x P x N

end