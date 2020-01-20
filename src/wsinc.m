function rf = wsinc(TBW, nt)
% wsinc -- Generate a Hamming windowed sinc pulse
%  Usage
%    rf = wsinc(TBW, nt)
%  Inputs
%    TBW           time-bandwidth product [sec] * [Hz = cycle/sec]
%    nt            number of samples
%  Outputs
%    rf            Hamming windowed sinc pulse
%

% Written by Namgyun Lee
% Email: ggang56@gmail.com
% Started: 04/13/2017, Last modified: 05/16/2018, 07/23/2018


% time-bandwidth product TBW = 2*N, total number of zeros
t = (-floor(nt/2):ceil(nt/2)-1).' ./ (floor(nt/2)); % time defined on [-1,1]

% when nt is odd, t is symmetric; otherwise, t has one more point in the
% negative side.
% nt = 4 => [-1, -0.5, 0, 0.5]
% nt = 3 => [-1, 0, 1]

% sinc(t)  : zero crossings at ..., -2  , -1  , 1  , 2  ,... defined on [-1,1]
% sinc(t*N): zero crossings at ..., -2/N, -1/N, 1/N, 2/N,... where N = TBW/2
% cos(pi*t): zero crossings at ..., -2  , -1  , 1  , 2  ,... defined on [-1,1]

% Calculate a Hamming windowed sinc pulse
rf = sinc(t * TBW / 2) .* (0.54 + 0.46 * cos(pi * t));

% Scale the waveform to sum to one
rf = rf / sum(rf);

end