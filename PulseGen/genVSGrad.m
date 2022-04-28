%% Generate gradient waveforms
%
% [G, t] = genVSGrad(T, pol, tend)
%
% in:
%      T    - struct of gradient and RF parameters
%      pol  - polarity of the VS gradients (e.g. [1,-1,1,-1])
%      tend - time at end of the gradient waveform
%
% out:
%      G       - gradient waveform with logical polarity
%      t       - gradient timecourse
%
% Written by Joseph G. Woods, CFMRI, UCSD, August 2020

function [G, t] = genVSGrad(T, pol, tend)

if ~exist('tend','var') || isempty(tend) % Default end time
    tend = T.RFe2; % To start of flip up RF
end

% If all polarity = 0, just generate an array of zeros
if ~any(pol)
    [G,t] = gengrad('t', 0, 0, 0, [0,tend], T.GUP);
    return;
end

% Ramp up raster formula
RUP_GRD_ms = @(A) round(ceil(round(round(A,12)*1e3/T.GUP,9))*T.GUP*1e-3, 3);

% Number of gradients to generate
nGrad = length(pol);

% Generate the gradient waveforms
gap = zeros(1, 2);
G   = zeros(round(tend*1e3/T.GUP), nGrad);
for ii = 1:nGrad
    gap(1)      = T.(['ta' num2str(ii)]);                      % Before G[ii]
    gap(2)      = RUP_GRD_ms( tend - T.(['td' num2str(ii)]) ); % Between G[ii] and end of flip up
    [G(:,ii),t] = gengrad('t', pol(ii)*T.Gmax, T.r, T.f, gap, T.GUP);
end
G = sum(G, 2); % Sum the gradient waveforms into one waveform

end