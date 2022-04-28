%% Generate gradient waveforms
% 
% [gout, tout] = gengrad(shape, Gmax, r, f, gap, dt)
%
% in:
%      shape - 'a'=linear ramp up, 'f'=flat top, 'd'=linear ramp down,
%              't'=trapezoid, 's'=sine, 'rect'=rectangle, 'c'=1-cosine
%      Gmax  - Maximum gradient amplitude (Units chosen by user: G/cm, mT/m, Hz/cm, etc)
%      r     - ramp time (ms)
%      f     - flat top time (ms)
%      gap   - Gaps between gradients (ms):
%              [before G1, after G1, after G2, after Gn-1, after Gn]
%      dt    - gradient raster time (µs)
%
% out:
%      gout - gradient waveform
%      tout - timepoints corresponding to gout (ms)
%
% Notes:
%      shape, Gmax, r, and f are arrays and should be the same length
%      gap is an array and should have one more element than shape
%
% Usage:
%      [gout, tout] = gengrad('tt', [1,-1], [r,r], [f,f], [0,1,0], GUP);
%
% Adapted from code by Jia Guo, UCSD, 2011
% Written by Joseph G. Woods, CFMRI, UCSD, May 2020

function [Gout,tout] = gengrad(shape, Gmax, r, f, gap, dt)

dt  = dt  * 1e-6; % convert temporal resolution from us to s
r   = r   * 1e-3; % convert ramp time from ms to s
f   = f   * 1e-3; % convert flat top time from ms to s
gap = gap * 1e-3; % convert gap from ms to s

nshape = numel(shape);  % number of gradient pulses
nr     = round(r/dt);   % ramp time resolution
nf     = round(f/dt);   % flat top time resolution
ngap   = round(gap/dt); % gap resoultion

% Need one more gap than shape
if nshape ~= numel(gap)-1
    error('It is required that nshape == numel(gap)-1.\n%s',...
         ['nshape = ' ns(nshape) ', nume(gap) = ' ns(ngap)])
end

% Insert first gap into gradient waveform
Gout = zeros(1, ngap(1));

% Loop over gradient pulses
for ii = 1:nshape
    
    % If shape(i) is empty, insert nothing
    G = [];
    
    switch shape(ii)
        
        case 'a' % linear ramp up
            G = Gmax(ii) * round( 1 : 1 : nr(ii) ) / nr(ii);
            
        case 'f' % flat top
            G = Gmax(ii) * ones(1, nf(ii));
            
        case 'd' % linear ramp down
            G = Gmax(ii) * round( nr(ii)-1 : -1 : 0 ) / nr(ii);
            
        case 't' % trapezoid
            ru = round( 1 : 1 : nr(ii) ) / nr(ii);
            ft = ones(1, nf(ii));
            rd = round( nr(ii)-1 : -1 : 0 ) / nr(ii);
            G  = Gmax(ii) * [ru, ft, rd];
            
        case 's' % sine
            ru = sin( round( 0 : 1 : floor(nf(ii)/2)  ) / floor(nf(ii)/2) * pi/2 );
            rd = ru(end:-1:1);
            G  = Gmax(ii) * [ru, rd];
            if length(G) < nf(ii)
                G = [G, zeros(1,nf(ii)-length(G))];
            end
            
        case 'r' % rectangle
            G = Gmax(ii) * ones(1, nf(ii));
            
        case 'c' % 1 - cosine
            ru = 1 - cos( round( 0 : 1 : floor(nf(ii)/2) ) / floor(nf(ii)/2) * pi/2 );
            rd = ru(end:-1:1);
            G  = Gmax(ii) * [ru, rd];
            if length(G) < nf(ii)
                G = [G, zeros(1,nf(ii)-length(G))];
            end
            
    end
    
    Gout = [Gout, G, zeros(1, ngap(ii+1))];
end

Gout = Gout(:);                                 % Gout as column vector
dt   = dt * 1e3;                                % dt in ms
tout = round(0 : dt : dt*(length(Gout)-1), 3)'; % tout in ms, round to nearest µs

end