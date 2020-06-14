function Mss = ssSSFP_BMF(theta, b1sqrd, TR, TE, tau, df, phi, F, lambda, tissuepars)
% ssSSFP_BMF -- Steady-state MT and flow bSSFP sequence using the extended matrix formalism
%  Usage
%     Mss = ssSSFP_BMF(theta, b1sqrd, TR, TE, tau, df, phi, F, lambda, tissuepars)
%  Inputs
%    theta        flip angle on resonance            [rad]
%    b1sqrd       mean square B1+ of a pulse         [uT^2]
%                 (over the duration of the pulse)
%    TR           repetition time                    [sec]
%    TE           echo time                          [sec]
%    tau          pulse duration                     [sec]
%    df           off-resonance                      [Hz]
%    phi          phase cycling increment            [rad]
%    F            blood flow                         [mL blood/g tissue/min]
%    lambda       tissue/blood partition coefficient [mL blood/g tissue]
%    tissuepars   structure containing all tissue parameters
%  Outputs
%    Mss          steady-state at TE (4 x 1)
%
% (c) Shaihan Malik 2019. King's College London
% Modified by Namgyun Lee, USC 09/04/2019


%% Unpack tissue parameters
M0f    = tissuepars.free.M0;  % equilibrium magnetiztion of the free pool
R1f    = tissuepars.free.R1;  % Free pool longitudinal relaxation rate [1/sec]
R2f    = tissuepars.free.R2;  % Free pool transverse relaxation rate   [1/sec]
M0s    = tissuepars.semi.M0;  % equilibrium magnetiztion of the semisolid pool
R1s    = tissuepars.semi.R1;  % Semisolid longitudinal relaxation rate [1/sec]
T2s    = tissuepars.semi.T2;  % semisolid T2, used in lineshape calculation [sec]
f      = tissuepars.f;        % the fractional size of the semisolid pool
k      = tissuepars.k;        % Overall exchange rate for free and both semisolid pools
G      = tissuepars.G;        % lineshape value at zero frequency [sec]
kf     = k * M0s;             % pseudo-first-order rate constant [1/sec]
ks     = k * M0f;             % pseudo-first-order rate constant [1/sec]

% gamma for RF calculation
% 1T = 1e4G = 1e6uT => 1G = 1e2uT
% [Hz/G] * [2*pi rad/cycle] * [G/1e2 uT] => [rad/sec/uT]
gam = 4257.746778 * 2 * pi * 1e-2; % [rad/sec/uT]

%% Calculate the mean saturation rate [rad/sec]
% [rad/sec/uT]^2 * [uT^2] * [sec] => [rad/sec]
W = pi * gam^2 * b1sqrd * G; % [rad/sec]

%% Calculate the inverse of Lambda + Gamma + Xi
Lambda = [0 0 0 0; 0 0 0 0; 0 0 -kf ks; 0 0 kf -ks]; % [1/sec]
Gamma  = [-F/lambda/60 0 0 0; 0 -F/lambda/60 0 0; 0 0 -F/lambda/60 0; 0 0 0 0]; % [1/sec]
Xi     = [-R2f 0 0 0; 0 -R2f 0 0; 0 0 -R1f 0; 0 0 0 -R1s]; % [1/sec]

invLGX = inv(Lambda + Gamma + Xi);

%% Calculate evolution operators
%--------------------------------------------------------------------------
% Calculate the exchange operator
%--------------------------------------------------------------------------
A = @(t) [f + 1     0                   0                            0                ;
            0     f + 1                 0                            0                ;
            0       0     1 + f * exp(-(f + 1) * ks * t)   1 - exp(-(f + 1) * ks * t) ;
            0       0     f - f * exp(-(f + 1) * ks * t)   f + exp(-(f + 1) * ks * t)] * 1 / (f + 1);

%--------------------------------------------------------------------------
% Clearance operator
% Describes the clearance of transverse and longitudinal magnetization by
% venous flow over a time interval [sec]
% F: [mL/g/min] * [min/60 sec] => [mL/g/sec]
%--------------------------------------------------------------------------
C = @(t) [exp(-F/60*t/lambda)              0                         0           0 ;
                    0              exp(-F/60*t/lambda)               0           0 ;
                    0                      0              exp(-F/60*t/lambda)    0 ;
                    0                      0                         0           1];

%--------------------------------------------------------------------------
% Calculate the relaxation operator
%--------------------------------------------------------------------------
E = @(t) diag([exp(-R2f * t); exp(-R2f * t); exp(-R1f * t); exp(-R1s * t)]);

%--------------------------------------------------------------------------
% Excitation is captured in a matrix Rx that contains a rotation part for
% the free pool and a saturation term for the restricted pool
% Rx: Rotation about the x-axis by theta
%--------------------------------------------------------------------------
R = [1         0             0             0      ;
     0    cos(theta)    sin(theta)         0      ;
     0   -sin(theta)    cos(theta)         0      ;
     0         0             0      exp(-W * tau)];

%--------------------------------------------------------------------------
% Calculate the free precession operator
% [rad] * [Hz] * [sec] => [rad]
%--------------------------------------------------------------------------
P = @(t) [ cos(2 * pi * df * t)    sin(2 * pi * df * t)    0    0 ;
          -sin(2 * pi * df * t)    cos(2 * pi * df * t)    0    0 ;
                      0                       0            1    0 ;
                      0                       0            0    1];

%--------------------------------------------------------------------------
% RF phase alternation
%--------------------------------------------------------------------------
Psi_bssfp = [ cos(phi)    sin(phi)    0    0 ;
             -sin(phi)    cos(phi)    0    0 ;
                 0           0        1    0 ;
                 0           0        0    1];

%--------------------------------------------------------------------------
% 4 x 4 identity matrix
%--------------------------------------------------------------------------
I = eye(4);

%--------------------------------------------------------------------------
% Define D
%--------------------------------------------------------------------------
D = [0; 0; (R1f + F / lambda /60) * M0f; R1s * M0s]; % [magnetiztion/g tissue/sec]

%% Now compute the steady-state solution
A1 = A(TR);
A2 = A(TE);

C1 = C(TR);
C2 = C(TE);

E1 = E(TR);
E2 = E(TE);

P1 = P(TR);
P2 = P(TE);

ACEP1 = A1 * C1 * E1 * P1;
ACEP2 = A2 * C2 * E2 * P2;

Mss = ACEP2 * inv(I - R * Psi_bssfp * ACEP1) * R * Psi_bssfp * (ACEP1 - I) * invLGX * D + (ACEP2 - I) * invLGX * D; % Mss(TE)

end