%% Returns the gyromagnetic ratio in either Hz/T or Hz/G.
%
% gam = gyroratio(type)
%
% in:
%      type -  'Hz/T', 'Hz/G', 'rad/s/T' or 'rad/s/G'
%
% out:
%      gam - gyromagnetic ratio
%
% Reference: https://physics.nist.gov/cgi-bin/cuu/Value?gammapbar
%
% Written by Joseph G. Woods, CFMRI, UCSD, April 2020

function gam = gyroratio(type)

gamMhzT = 42.577478518; % In MHz/T

switch type
    
    case 'Hz/T'
        gam = gamMhzT * 1e6;
    
    case 'Hz/G'
        gam = gamMhzT * 1e2;
        
    case 'rad/s/T'
        gam = gamMhzT * 2*pi * 1e6;
    
    case 'rad/s/G'
        gam = gamMhzT * 2*pi * 1e2;
        
    otherwise
        error('type should be one of: Hz/T, Hz/G, rad/s/T, rad/s/G')
        
end
