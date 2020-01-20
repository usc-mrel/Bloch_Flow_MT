function m = bloch_flow(RFpulses, TR, TE, T1, T2, df, f, lambda, m0, asl_bolus)
% BLOCH_FLOW -- Bloch equations with flow (perfusion)
%  Usage
%    m = bloch_flow(RFpulses, TR, TE, T1, T2, df, f, lambda, m0, asl_bolus)
%  Inputs
%    RFpulses    complex RF pulse train (nt x 1)
%    TR          TR [sec] (nt x 1)
%    TE          TE [sec] (nt x 1)
%    T1          T1 [sec] (1 x 1)
%    T2          T2 [sec] (1 x 1)
%    df          tissue off-resonance frequency [Hz] (1 x 1)
%    f           blood flow [mL/g/min] (1 x 1)
%    lambda      muscle/blood partition coefficient for water [mL of blood /g of myocardium]
%    m0          equilibrium z magnetization
%    asl_bolus   ASL bolus signal for each interval (nt x 1)
%  Outputs
%    m           magnetization vector (3 x nt)
%

%% Define function handles
I = eye(3); % 3 x 3 identity matrix

% Rotations are left-handed
% Rotation by alpha about the x-axis
Rx = @(alpha) [     1            0            0      ;
                    0       cos(alpha)    sin(alpha) ;
                    0      -sin(alpha)    cos(alpha)];

Ry = @(alpha) [ cos(alpha)      0        -sin(alpha) ;
                    0           1              0     ;
                sin(alpha)      0         cos(alpha)];

Rz = @(alpha) [ cos(alpha)   sin(alpha)        0     ;
               -sin(alpha)   cos(alpha)        0     ;
                    0           0              1    ];

% free precession over a period tau about the z-axis
% [2*pi rad/cycle] * [Hz] * [sec]
P = @(tau) [ cos(2*pi*df*tau)  sin(2*pi*df*tau)         0 ;
            -sin(2*pi*df*tau)  cos(2*pi*df*tau)         0 ;
                    0                 0                 1];

% T1 and T2 relaxation over a time tau
E = @(tau) [exp(-tau/T2)        0              0     ;
                  0       exp(-tau/T2)         0     ;
                  0             0       exp(-tau/T1)];

% Clearance by venous flow: blood flow [mL/g/min] * [min/60 sec] => [mL/g/sec]
C = @(tau) [exp(-f/60*tau/lambda)            0                      0           ;
                     0             exp(-f/60*tau/lambda)            0           ;
                     0                       0            exp(-f/60*tau/lambda)];

% T1 Recovery + clearance by venous flow
D = @(tau) (I - C(tau) * E(tau)) * [0 0 m0].';

%% Perform Bloch simulations with flow
nt = length(RFpulses); % number of time intervals
m = zeros(3,nt, 'double');
ma = [0 0 m0].';

alpha = abs(RFpulses);   % RF pulse angle [rad]
phi   = angle(RFpulses); % RF pulse phase [rad]

for k = 1:nt
    %======================================================================
    %      RF(1)         RF(2)         RF(nt-1)       RF(nt)
    %        |             |             |             |
    %     ma | mb mc    md |             |             |
    %       \|/   |       \|             |             |
    %   -----+----x--------+----x--------+----x--------+----x--------
    %        t(1)          t(2)          t(nt-1)       t(nt)
    %        <------------><------------><------------><------------>
    %           TR(1)         TR(2)         TR(nt-1)      TR(nt)
    %        <--->         <--->         <--->         <--->
    %        TE(1)         TE(2)         TE(nt-1)      TE(nt)
    %======================================================================
    TR1 = TR(k);
    TE1 = TE(k);
    alpha1 = alpha(k); % RF pulse angle [rad]
    phi1   = phi(k);   % RF pulse phase [rad]

    % Calculate the amount of arterial blood flowing into the voxel during
    % the kth interval using the hard pulse approximation
    bolus1 = asl_bolus(k);

    % P: free precession over a period tau about the z-axis
    % E: T1 and T2 relaxation over a time tau
    % C: clearance by venous flow
    % D: T1 Recovery + clearance by venous flow
    ma = ma + [0 0 bolus1 * TR1].'; % inflow of the tagged blood
    mb = Rz(-phi1) * Rx(alpha1) * Rz(phi1) * ma;
    mc = C(TE1) * P(TE1) * E(TE1) * mb + D(TE1);
    md = C(TR1 - TE1) * P(TR1 - TE1) * E(TR1 - TE1) * mc + D(TR1 - TE1);
    m(:,k) = mc;
    ma = md;
end

end