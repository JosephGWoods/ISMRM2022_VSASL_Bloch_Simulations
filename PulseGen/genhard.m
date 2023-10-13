%% A function to generate a hard RF pulse
%
% [B1, t_ms, bw] = genhard(FA, phase, dur, B1max, dt, units)
%
% in:
%      FA    - flip angle of B1+ pulse (degrees)
%      phase - phase of B1+ pulse (degrees)
%      dur   - duration of B1+ pulse (ms)
%      B1max - amplitude of B1+ pulse (match with units)
%      dt    - RF pulse step size or resolution (µs)
%      units - 'G' Gauss, 'T' Tesla, 'Hz' Hertz
%
% out:
%      B1    - complex B1+ waveform (Gauss - default)
%      t     - time array (ms)
%      bw    - bandwidth of RF pulse (Hz)
%
% NB:  If both dur and B1max are supplied, the duration will be used and
%      the appropriate amplitude will be calculated. If this amplitude is
%      larger than the supplied B1max, then an error is thrown.
%
% Rectangular pulse example:
%
%      B1+
%    ___|___
%   |   |   |     
% __|___|___|__ t
%    <----->
%       T
%
% B1+ = FA / γT  ,  |t| ≤ T/2
%
% Written by Joseph G. Woods, CFMRI, UCSD, May 2020

function [B1, t, bw, dur] = genhard(FA, phase, dur, B1max, dt, units)

if nargin < 6
    error('Not enough input arguments!')
end

switch units
    case 'G' ; gamrad = gyroratio('rad/s/G');
    case 'T' ; gamrad = gyroratio('rad/s/T');
    case 'Hz'; gamrad = 2*pi; % Simply do not divide by γ in B1max
    otherwise; error('Units can only be G, T, or Hz!');
end

FA     = deg2rad(FA   ); % Convert to radians
phase  = deg2rad(phase); % Convert to radians

% Check if pulse duration and amplitude are supplied
if isempty(dur) && ~isempty(B1max)
    dur = 1e-3 * ceil( 1e6 * FA / (gamrad * B1max) / dt ) * dt; % In ms
elseif ~isempty(dur) && isempty(B1max)
    B1max = inf;
elseif isempty(dur) && isempty(B1max)
    error('User must supply either pulse duration or amplitude!')
end

dur_us = dur * 1e3;        % Convert to µs
N      = round(dur_us/dt); % Number of support points
t_us   = 0:dt:dt*(N-1);    % Waveform time
t      = t_us * 1e-3;      % Convert to ms for output only
bw     = 1/(dur*1e-3);     % Calculate the bandwidth of the RF pulse in Hz

%% Calculate B1+ waveform

% Form rectangular pulse envelope
B1 = ones(N, 1);

%% Calculate B1max required to achieve flip angle

% Pulse integral
area = dt * 1e-6 * sum(B1);

% B1max calculation (hard pulse FA scales linearly with B1max)
B1maxCalc = FA / ( gamrad * area );
if B1maxCalc > B1max; error('B1 of hard pulse is > B1max'); end

% Scale RF pulse to achieve flip angle
B1 = B1 * B1maxCalc;

%% Add constant phase to B1+ pulse

B1 = B1 .* exp(1i * phase);

