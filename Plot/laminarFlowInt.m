%% Calculate laminar flow integration of Mz for a given set of velocities
%
% in:
%      Mz - array of Mz values at each velocity (arbitrary units)
%      v  - array of velocities (arbitrary units)
%
% out:
%      Mz_laminar - Mz integrated over velocities
%      lind       - indicies of mean velocities
%
% Written by Joseph G. Woods, University of Oxford, April 2022
% Edited by Dapeng Liu, Johns Hopkins University, May 2022

function [Mz_laminar, lind] = laminarFlowInt(Mz, v)

% Find zero velocity
zvind = find(v==0);

nv1 = zvind-findNearest(min(v)/2,v); % number of -ve mean velocities we can do the laminar integration for [min(v)/2]
nv2 = findNearest(max(v)/2,v)-zvind; % number of +ve mean velocities we can do the laminar integration for [max(v)/2]

% Laminar flow integration (uniform distribution from 0 to Vmax)
Mz_laminar = zeros(1,nv1+nv2-1);

% Negative velocities
for ii = 0:nv1-1
    ind = findNearest(2*v(zvind-ii),v);         % find max velocity for current mean velocity (2*meanV)
    Mz_laminar(nv1-ii) = mean(Mz(zvind:-1:ind)); % mean Mz across laminar flow profile
end

% Positive velocities
for ii = 1:nv2-1
    ind = findNearest(2*v(zvind+ii),v);        % find max velocity for current mean velocity (2*meanV)
    Mz_laminar(nv1+ii) = mean(Mz(zvind:1:ind)); % mean Mz across laminar flow profile
end

% Indices of mean velocities
lind = zvind-nv1+1:zvind+nv2-1;

end