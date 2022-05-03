%% A function to generate a windowed sinc RF pulse
%
% B1 = gensinc(FA, phase, dur, NL, NR, alpha, dt, units, slthk)
%
% in:
%      FA     - flip angle of B1+ pulse (degrees)
%      phase  - phase of B1+ pulse (degrees)
%      dur    - duration of B1+ pulse (ms)
%      NL     - number of zero crossings on left of pulse
%      NR     - number of zero crossings on right of pulse
%      alpha  - cosine window tuning parameter
%      dt     - RF pulse step size or resolution (µs)
%      units  - 'G' Gauss (default), 'T' Tesla, 'Hz' Hertz
%
% out:
%      B1    - complex B1+ waveform (units)
%
% notes:
%        TBWP:  NL + NR is the total number of zero crossings (time-bandwidth
%               product)
%        alpha: α = 0 is truncated sinc (no windowing)
%               α = 0.5 is Hanning window (continuous 1st derivative for
%               NL = NR)
%               α = 0.46 is Hamming window (12.5 reduced 1st derivative at
%               margins for NL = NR)
%
%
% Example: sinc pulse with TBWP = 4
%
% B1+     __
% |      /  \
% |_____/____\____ t
% |  \_/      \_/
% |  <-><-><-><->
%    t0 t0  t0 t0
%
% B1(t) = { B1max*t0 [(1-α)+αcos(πt/Nt0)] sin(πt/t0)/πt , -NLt0 ≤ t ≤ NRt0
%         { 0                                           ,  elsewhere
%
% Reference: Bernstein et al., Handbook of MRI Pulse Sequences (2004) page 39
%
% Written by Joseph G. Woods, CFMRI, UCSD, May 2020

function B1 = gensinc(FA, phase, dur, NL, NR, alpha, dt, units)

if nargin < 3
    error('Not enough input arguments! Must supply FA, phase, and dur!')
end
if ~exist('NL'    ,'var')||isempty(NL    ); NL     =  2.0     ; end
if ~exist('NR'    ,'var')||isempty(NR    ); NR     =  2.0     ; end
if ~exist('alpha' ,'var')||isempty(alpha ); alpha  = 'Hamming'; end
if ~exist('dt'    ,'var')||isempty(dt    ); dt     =  2.0     ; end
if ~exist('units' ,'var')||isempty(units ); units  =  'G'     ; end

switch alpha
    case 'Hann';    alpha = 0.5;
    case 'Hamming'; alpha = 0.46;
end

switch units
    case 'G' ; gam = gyroratio('Hz/G');
    case 'T' ; gam = gyroratio('Hz/T');
    case 'Hz'; gam = 1; % Simply do not divide by γ in B1max
end

FA     = deg2rad(FA   );        % Convert to radians
phase  = deg2rad(phase);        % Convert to radians
dur_us = dur * 1e3;             % Convert to µs
N      = max( NL, NR );         % N is the larger of NL and NR
t0     = dur_us / (N*2);        % t0 is duration of each lobe (centre lobe is 2*t0)
t_us   = ((-NL*t0+dt/2):dt:(NR*t0-dt/2))'; % Waveform time

%% Calculate B1+ waveform

% Sinc B1+ envelope formula
B1 = t0 * ( (1-alpha) + alpha*cos(pi*t_us/(N*t0)) ) .* sin(pi*t_us/t0) ./ (pi*t_us);

% B1+ is nan at t=0 (centre of the pulse)
B1(t_us==0) = 1;

%% Calculate B1max required to achieve flip angle

% Pulse integral
area = dt * 1e-6 * sum( sqrt(real(B1).^2 + imag(B1).^2) ); % The imag part is not needed since B1 is real-valued here.

% B1max calculation (sinc pulse FA scales linearly with B1max)
B1max = FA / ( 2 * pi * gam * area );

% Scale sinc pulse to achieve flip angle
B1 = B1 * B1max;

%% Add constant phase

% Add phase to B1+ pulse
B1 = B1 .* exp(1i * phase);

end
