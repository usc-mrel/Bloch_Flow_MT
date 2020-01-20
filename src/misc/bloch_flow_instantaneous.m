function m = bloch_flow_instantaneous(RFpulses, TRs, TEs, T1, T2, df, F, lambda, M0, ASL_bolus)
% bloch_flow_instantaneous -- Bloch-Flow simulator using the extended matrix formalism
%                             assuming instantaneous RF excitation
%  Usage
%    m = bloch_flow_instantaneous(RFpulses, TRs, TEs, T1, T2, df, F, lambda, m0, ASL_bolus)
%  Inputs
%    RFpulses    vector of complex RF pulses        [rad]      (N x 1)
%    TRs         vector of repetition times         [sec]      (N x 1)
%    TEs         vector of echo times               [sec]      (N x 1)
%    T1          longitudinal relaxation time       [sec]      (1 x 1)
%    T2          transverse relaxation time         [sec]      (1 x 1)
%    df          off-resonance                      [Hz]       (1 x 1)
%    F           blood flow                         [mL/g/min] (1 x 1)
%    lambda      tissue/blood partition coefficient [mL blood/g tissue]          (1 x 1)
%    M0          equilibrium magnetization          [magnetization/g tissue]     (1 x 1)
%    ASL_bolus   ASL bolus signal                   [magnetization/g tissue/sec] (N x 1)
%  Outputs
%    m           magnetization vector               [magnetization/g tissue] (3 x N)
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
D = @(tau) (I - C(tau) * E(tau)) * [0 0 M0].';

%% Perform Bloch-Flow simulation
N = length(RFpulses); % number of time intervals
m = zeros(3,N, 'double');
ma = [0 0 M0].';

alpha = abs(RFpulses);   % RF pulse angle [rad]
phi   = angle(RFpulses); % RF pulse phase [rad]

for k = 1:N
    %======================================================================
    % Timing diagram
    %----------------------------------------------------------------------
    %      RF(1)         RF(2)         RF(N-1)       RF(N)
    %        |             |             |             |
    %     ma | mb mc    md |             |             |
    %       \|/   |       \|             |             |
    %   -----+----x--------+----x--------+----x--------+----x--------
    %        t(1)          t(2)          t(N-1)        t(N)
    %        <------------><------------><------------><------------>
    %           TR(1)         TR(2)         TR(N-1)        TR(N)
    %        <--->         <--->         <--->         <--->
    %        TE(1)         TE(2)         TE(N-1)       TE(N)
    %======================================================================
    TR = TRs(k); % current repetition time
    TE = TEs(k); % current echo time

    %----------------------------------------------------------------------
    % Labeled flow
    % Calculate the amount of labeled inflow during the kth interval using 
    % the hard pulse approximation
    %----------------------------------------------------------------------
    ma = ma + [0 0 ASL_bolus(k) * TR].'; % inflow of labeled blood

    %----------------------------------------------------------------------
    % RF excitation
    %----------------------------------------------------------------------
    mb = Rz(-phi(k)) * Rx(alpha(k)) * Rz(phi(k)) * ma;

    %----------------------------------------------------------------------
    % Relaxation (E) - Free precession (P) - Flow (F)
    %----------------------------------------------------------------------
    mc = C(TE) * P(TE) * E(TE) * mb + D(TE);

    %----------------------------------------------------------------------
    % Relaxation (E) - Free precession (P) - Flow (F)
    %----------------------------------------------------------------------
    md = C(TR - TE) * P(TR - TE) * E(TR - TE) * mc + D(TR - TE);
    m(:,k) = mc;
    ma = md;
end

end