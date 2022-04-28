% A function to finds nearest number to the one specified in the array
% provided.
%
% [ind,res] = findNearest(x, y, n, d)
%
% in:
%      x - value to search for in y (scaler, xL-by-1, or 1-by-xL)
%      y - array within which to search for n (array)
%      n - output n nearest elements to each x (scalar)
%      d - maximum difference between n and values in y
%
% out:
%      ind - linear indice(s) in y for the nearest element to n (xL by n)
%      res - residual difference between nearest y element and n (xL by n)
%
% Written by Joseph G. Woods, FMRIB, Oxford, March 2019

function [ind,res] = findNearest(x, y, n, d)

% Number of nearest elements is 1 by default
if ~exist('n','var') || isempty(n)
    n = 1;
end

% Maximum difference is inf by default
if ~exist('d','var') || isempty(d)
    d = inf;
end

% Initialise outputs (scalars if scalar n)
res = zeros(n, length(x));
ind = zeros(n, length(x));

for ii = 1:1:length(x)
    
    minNum       = abs( y - x(ii) );       % Calculate distances
    [rtmp, itmp] = sort(minNum, 'ascend'); % Rank distances from smallest to largest
    ind(:,ii)    = itmp(1:n);              % Take the top n elements
    res(:,ii)    = rtmp(1:n);              % Take the top n elements
    
end

if max(res) > d
    warning('The max(nearest value) is %d from the requested number. Consider investigating', max(res));
end

end

