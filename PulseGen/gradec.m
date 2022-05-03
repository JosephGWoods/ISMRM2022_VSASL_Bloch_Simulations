%% Eddy current induced gradient fields
%
% gec = gradec(G, A, tau, GUP)
%
% in:
%      G     - Input gradient waveform (any units)
%      A     - Array of amplitude of eddy-currents (fraction)
%      tau   - Array of time-constants of eddy-currents (ms)
%      GUP   - Gradient update time (µs)
%
% out:
%      gec  - Modelled additional gradient field due to eddy-currents (same units as G)
%
% Notes:
%      G is an array with size Nt x 1
%      A has size 1 x Nτ
%      tau has size 1 x Nτ
%      GUP is a scalar
%      gec has size Nt x Nτ
%
% This function numerically calculates the additional gradient field
% induced by switching of the gradient amplitude. This is modelled as a sum
% of exponentials. See Vaals and Bergman, JMR 1990
% (https://doi.org/10.1016/0022-2364(90)90365-G)
%
% Equation:
%           gec(t,n) = -dG/dt * H(t) A exp( -t / τ )
%
% where * is a convolution and H is the unit step-function (heaviside).
% The gradient is calculated on a per gradient raster interval and uses
% finite differences.
%
% Written by Joseph G. Woods, CFMRI, UCSD, May 2020

function Gec = gradec(G, A, tau, GUP)

tau = tau * 1e3; % Convert from ms to µs

% Save size of inputs
sizg   = size(G);
sizA   = size(A);
siztau = size(tau);

% Convert inputs to matrices
G   = reshape(G  , sizg(1), []); tL   = sizg(1); GL = size(G,2);
A   = reshape(A  , 1      , []); AL   = length(A);
tau = reshape(tau, 1      , []); tauL = length(tau);

% Check dimensions of inputs
if AL>1 && tauL>1 && AL~=tauL
    error('The size of A and tau must be the same: size(A) = %s, size(tau) = %s.\n',...
        array2string(sizA),array2string(siztau));
end
if GL>1 && AL>1 && GL~=AL
    error('The non-time size of g and A must be the same: size(g,2:end) = %s, size(A) = %s.\n',...
        array2string(sizg(2:end)),array2string(sizA));
end
if GL>1 && tauL>1 && GL~=tauL
    error('The non-time size of g and tau must be the same: size(g,2:end) = %s, size(tau) = %s.\n',...
        array2string(sizg(2:end)),array2string(siztau));
end
if length(GUP)>1; error('GUP should be a scalar'); end

% Create the time array for this gradient
t = (GUP : GUP : GUP*tL)'; % (µs)

% Repmat vectors to have the same size
if GL  >1; t   = repmat(t  , 1 , GL); end
if AL  >1; A   = repmat(A  , tL, 1 ); end
if tauL>1; tau = repmat(tau, tL, 1 ); end

% Calculate the rate of change of the gradient waveform per gradient raster
% point.
G    = cat(1, zeros(1, GL), G); % Pad g with leading zero
dGdt = diff(G, [], 1);          % Finite difference through time

% Most time-consuming part of code
% ec  = heaviside(t) .* A .* exp(-t./tau); % Leave out heaviside(t) since this simply = 1 for all t>0 which everything is by default
ec  = A .* exp(-t./tau);                   % Calculate the eddy-current modulation function
Gec = zeros(tL*2-1, max(AL,tauL));         % Preallocate for speed
for ii = 1:max(AL,tauL)                    % Faster to not use parfor due to overheads % arrayfun slightly slower than for loop
    Gec(:,ii) = conv(-dGdt, ec(:,ii));     % Calculate the additional gradient due to these eddy-currents
end

% Alternative: convmtx should be faster but is 10x slower
% ec  = A .* exp(-t./tau);
% c   = convmtx(-dgdt, tL);
% gec = c * ec;

% Alternative: ifft(fft(A).*fft(B)) can be faster than conv(A,B), but not here
% ec   = A .* exp(-t./tau);
% dgdt = [dgdt; zeros(tL-1,1)];
% ec   = [ec; zeros(tL-1,max(AL,tauL))];
% gec  = zeros(tL*2-1, max(AL,tauL));
% for ii = 1:max(AL,tauL)
%     gec(:,ii) = ifft(fft(-dgdt).*fft(ec(:,ii)));
% end

% Select only the first half of the convolution
Gec = Gec(1:tL,:);

end