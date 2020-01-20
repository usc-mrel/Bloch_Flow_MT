function [mx,my,mz] = bloch_flow_MA(b1, taus, zetas, T1, T2, df, mode, mx0, my0, mz0, F, lambda, m0, ASL_bolus)
% bloch_flow -- Bloch-Flow simulator without gradient using the extended matrix formalism
%               (similar to the minimal approximation MT style approach by Portnoy and Stanisz MRM 2007)
%  Usage
%    [mx,my,mz] = bloch_flow_MA(b1, taus, zetas, T1, T2, df, dp, mode, mx0, my0, mz0, F, lambda, m0, ASL_bolus)
%  Inputs
%    b1          RF pulse in G. Can be complex          [G]        (M x 1)
%    taus        time duration of each b1 and gr point  [sec]      (M x 1)
%    zetas       measurement time of each interval      [sec]      (M x 1)
%    T1          longitudinal relaxation time           [sec]      (1 x 1)
%    T2          transverse relaxation time             [sec]      (1 x 1)
%    df          off-resonance frequency                [Hz]       (1 x 1)
%    mode        Bitmask mode:
%                Bit 0: 0-Simulate from start or M0, 1-Steady State (not supported)
%                Bit 1: 1-Record m at time points. 0-just end time.
%    mx0         starting x magnetization                          (1 x 1)
%    my0         starting y magnetization                          (1 x 1)
%    mz0         starting z magnetization                          (1 x 1)
%    F           blood flow                             [mL/g/min] (1 x 1)
%    lambda      tissue/blood partition coefficient     [mL blood/g tissue]          (1 x 1)
%    m0          equilibrium magnetization              [magnetization/g tissue]     (1 x 1)
%    ASL_bolus   ASL bolus signal                       [magnetization/g tissue/sec] (M x 1)
%  Outputs
%    mx          array of the resulting x magnetization [magnetization/g tissue] (M x 1) or (1 x 1)
%    my          array of the resulting y magnetization [magnetization/g tissue] (M x 1) or (1 x 1)
%    mz          array of the resulting z magnetization [magnetization/g tissue] (M x 1) or (1 x 1)
%
% Copyright (c) 2019. Namgyun Lee


%% Define function handles
I = eye(3); % 3 x 3 identity matrix

%--------------------------------------------------------------------------
% Rotation matrices
%--------------------------------------------------------------------------
% Rotation by alpha about the x-axis
Rx = @(alpha) [     1              0             0     ;
                    0         cos(alpha)    sin(alpha) ;
                    0        -sin(alpha)    cos(alpha)];

% Rotation by alpha about the y-axis
Ry = @(alpha) [ cos(alpha)        0        -sin(alpha) ;
                    0             1              0     ;
                sin(alpha)        0         cos(alpha)];

% Rotation by alpha about the z-axis
Rz = @(alpha) [ cos(alpha)    sin(alpha)         0     ;
               -sin(alpha)    cos(alpha)         0     ;
                    0            0               1    ];

%--------------------------------------------------------------------------
% free precession over a period tau about the z-axis
%--------------------------------------------------------------------------
% [2*pi rad/cycle] * [Hz] * [sec] => [rad]
P = @(tau) [ cos(2*pi*df*tau)    sin(2*pi*df*tau)            0 ;
            -sin(2*pi*df*tau)    cos(2*pi*df*tau)            0 ;
                     0                   0                   1];

%--------------------------------------------------------------------------
% Relaxation operator
% Describes T1 and T2 relaxation over a time interval of duration "tau" [sec]
%--------------------------------------------------------------------------
E = @(tau) [exp(-tau/T2)          0               0      ;
                  0         exp(-tau/T2)          0      ;
                  0               0         exp(-tau/T1)];

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

m_zeta = zeros(M,3, 'double');

%--------------------------------------------------------------------------
% Set up the starting magnetization
%--------------------------------------------------------------------------
ma = [mx0 my0 mz0].';

for idx1 = 1:M
    %=======================================================================
    % Timing diagram
    %----------------------------------------------------------------------
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
    %=======================================================================
    tau  = taus(idx1);  % current duration of the interval
    zeta = zetas(idx1); % current measurement time

    %----------------------------------------------------------------------
    % RF excitation
    %----------------------------------------------------------------------
    b1x = real(b1(idx1)); % [G]
    b1y = imag(b1(idx1)); % [G]
    theta = gam * 1e2 * sqrt(b1x * b1x + b1y * b1y) * tau; % [rad/sec/uT] * [1e2uT/G] * [G] * [sec] => [rad]
    phi = atan2(b1y, b1x); % [rad]
    mb = Rz(-phi) * Rx(theta) * Rz(phi) * ma;

    %----------------------------------------------------------------------
    % Labeled flow
    % Calculate a new pseudo M0 term: M0 + s(t) * T1app, which is the 
    % hypothetical equilibrium magnetization for modeling flow during this
    % short interval.
    %----------------------------------------------------------------------
    pseudo_M0 = [0 0 m0 + ASL_bolus(idx1) * T1app].'; % 3 x 1

    %----------------------------------------------------------------------
    % Relaxation (E) - Precession (P) - Flow (F)
    %----------------------------------------------------------------------
    CPE = C(zeta) * P(zeta) * E(zeta);
    mc = CPE * mb + (I - CPE) * pseudo_M0;

    %----------------------------------------------------------------------
    % Record the magnetization at zeta
    %----------------------------------------------------------------------
    m_zeta(idx1,:) = mc;

    %----------------------------------------------------------------------
    % Relaxation (E) - Precession (P) - Flow (F)
    %----------------------------------------------------------------------
    CPE = C(tau - zeta) * P(tau - zeta) * E(tau - zeta);
    md = CPE * mc + (I - CPE) * pseudo_M0;
    ma = md;
end

if mode == 2 % Record m at time points
    mx = m_zeta(:,1); % M x 1
    my = m_zeta(:,2); % M x 1
    mz = m_zeta(:,3); % M x 1
elseif mode == 0 % just end time
    mx = m_zeta(M,1); % 1 x 1
    my = m_zeta(M,2); % 1 x 1
    mz = m_zeta(M,3); % 1 x 1
end

end