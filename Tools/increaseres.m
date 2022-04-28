% A function to increase the resolution of a 1-D array by repeating
% elements contiguously.
%
% y = increaseres(x, m, dim)
%
% in:
%      x   - input 1-D array
%      m   - integer factor to increase resolution by
%      dim - dimension to increase resolution along
%
% out:
%      y   - Increased resolution array
%
% Written by Joseph G. Woods, CFMRI, UCSD, April 2020

function y = increaseres(x, m, dim)

if ndp(m)~=0
    error('m must be an integer')
end

% Find the first non-singular dimension
if ~exist('dim','var') || isempty(dim)
    dim = fnsd(x);
end

% Intialise y with the new array size
sizx      = size(x);
sizy      = sizx;
sizy(dim) = sizy(dim) * m;
y         = zeros(sizy);

% Increase the resolution of x
for ii = 1:m
    y(ii:m:end) = x; % FIXME: only works for 1 dimension presently
end

end
