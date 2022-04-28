%% Find the first non-singular dimension
%
% dim = fnsd(x)
%
% in:
%      x   - input n-D array
%
% out:
%      dim - The first non-singular dimension
%
% Written by Joseph G. Woods, CFMRI, UCSD, May 2020

function dim = fnsd(x)

dim = find(size(x)~=1);

if ~isempty(dim)
    dim = dim(1);
else
    dim = 1;
end

end