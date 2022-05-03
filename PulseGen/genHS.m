% Non-selective hyperbolic-secant pulse morphed into a cosine envelope ala VERSE
%
% [rho, theta, res] = genHS(beta, mu, dtus, cexp)
%
% in:
%      beta   - modulation angular frequency (rad/s)
%      mu     - shaping parameter (slab width to tranisition width ratio)
%      dtus   - time resolution (µs)
%      cexp   - output envelope is cos^cexp
%      initph - initial phase (rad)
%
% out:
%      rho   - RF envelope (normalised)
%      theta - RF phase (in rad)
%      res   - RF resolution
%
% Originally written by Eric C. Wong CFMRI, UCSD
% Edited by Joseph G. Woods, CFMRI, UCSD, May 2020

function [rho, theta, res] = genHS(beta, mu, dtus, cexp, initph)

dt     = 1e-6 * dtus;         % convert to seconds
initph = initph * pi / 180.0; % convert to radians

% area of cos^cexp from 0 to pi/2
ac = (sqrt(pi)/cexp) * exp(gammaln(0.5*(1.+cexp))-gammaln(0.5*cexp));

% this duration makes the area of the pulse match that of the sech
dur = pi*pi/(4*beta*ac);
n   = floor(dur/dt);
dp  = pi/(2*n);
res = 2*n;
cint_scale = 0.5*pi*dp/ac;

% Calculate waveforms - right side only, then replicate left
% start cos^cexp integral by integrating half way into the first time bin
cint = 0.5*cint_scale;

rho   = zeros(res, 1);
theta = zeros(res, 1);
for ii = 0 : 1 : n-1
    
    % envelope is cos^cexp
    p = (ii+0.5)*dp;
    rho(n+ii+1) = cos(p)^cexp;
    rho(n-ii)   = rho(n+ii+1);
    
    % inverse of integral of sech is 2*atanh(tan(x/2))
    bt            =  2 * atanh( tan(0.5*cint) );
    sech          =  2 / (exp(bt) + exp(-1*bt));
    phase         = mu * log(sech) + initph;
    theta(n+ii+1) = modp(phase);
    theta(n-ii)   = theta(n+ii+1);
    
    % integrate cos^cexp and normalize to [0,pi/2] to find bt to use in phase calc
    cint = cint + cint_scale * cos( (ii+1)*dp )^cexp;
    
end

    % Wraps phase into -pi to pi
    function b = modp(a)
        if a < 0; b = (2 * rem((a-pi)/2,pi) + pi);
        else;     b = (2 * rem((a+pi)/2,pi) - pi);
        end
    end

end
