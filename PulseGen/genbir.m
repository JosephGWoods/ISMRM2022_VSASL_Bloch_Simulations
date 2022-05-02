% BIR segments
%
% [rho, theta, nbir] = genbir(wmax, zeta, tkap, pw, dtus)
%
% in:
%      wmax - max frequency sweep in bir module (hz)
%      zeta - zeta in bir module (s^-1)
%      tkap - tan of kappa in bir module
%      pw   - pulse width of adiabatic half passage (ms)
%             (full passage is twice as long)
%      dtus - time resolution (µs)
%
% out:
%      rho   - struct of RF envelope segments (normalised)
%      theta - struct of RF phase segments (in rad)
%      nbir  - RF resolution
%
% Originally written by Eric C. Wong CFMRI, UCSD
% Edited by Joseph G. Woods, CFMRI, UCSD, June 2020

function [rho, theta, nbir] = genbir(wmax, zeta, tkap, pw, dtus)

dt   = 1.e-6 * dtus;         % convert to seconds
dur  = 1.e-3 * pw;           % convert to seconds
kap  = atan(tkap);
w    = wmax * 2 * pi * dt;   % convert to phase per dt
n    = round(dur/dt);

% Initialise the arrays
rho.birleft    = zeros(  n, 1);
rho.birmid     = zeros(2*n, 1);
rho.birright   = zeros(  n, 1);
theta.birleft  = zeros(  n, 1);
theta.birmid   = zeros(2*n, 1);
theta.birright = zeros(  n, 1);

% Generate RF amplitude and phase
ptmp = 0.0;
for ii = 0:(n-1)
    
    t = ii / (n-1);
    
    rho.birleft(ii+1)   = tanh(zeta*(1-t)); % RF amplitude
    theta.birleft(ii+1) = modp(ptmp);       % RF phase
    theta.birmid(n-ii)  = modp(ptmp+pi);    % RF phase
    
    pstep = w * tan(kap*t) / tkap;
    ptmp  = ptmp + pstep;
    
end
nbir = n;


for ii = 0:(n-1)
    rho.birmid(ii+1)     = rho.birleft(n-ii);
    rho.birmid(n+ii+1)   = rho.birleft(ii+1);
    rho.birright(ii+1)   = rho.birleft(n-ii);
    theta.birmid(n+ii+1) = theta.birmid(n-ii);
    theta.birright(ii+1) = theta.birleft(n-ii);
end

    % Wraps phase into -pi to pi
    function b = modp(a)
        if a < 0; b = (2 * rem((a-pi)/2,pi) + pi);
        else;     b = (2 * rem((a+pi)/2,pi) - pi);
        end
    end

end
