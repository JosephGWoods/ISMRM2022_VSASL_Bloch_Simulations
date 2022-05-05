%% Function to generate a BIR-4 velocity selective module
%
% [B1, GLabel, GCont, T] = genBIR4(T, bSection)
%
% in:
%      T            - struct of gradient and RF parameters 
%      bSection     - sections of the module to generate: 'excite',
%                     'refocus', 'VSgrad', 'combine', or 'all'
%      bvelCompCont - flag to use velocity-compensated control
%
% out:
%      B1     - complex B1+ waveform (units of B1max)
%      GLabel - label module gradient waveform  (units of Gmax)
%      GCont  - control module gradient waveform (units of Gmax)
%      T      - updated struct of gradient and RF parameters
%
% T parameter descriptions:
%      B1max  - B1+ maximum amplitdue (units)
%      Gmax   - maximum gradient amplitude (units/cm)
%      RFUP   - RF update time (µs)
%      GUP    - gradient update time (µs)
%      units  - B1+ and gradient units ('G', 'T', 'Hz')
%      f     - gradient flat top time
%      r     - gradient rise time
%      ta[n] - start of trapezoid attack time for nth gradient
%      td[n] - end of trapezoid decay time for nth gradient
%      RFe   - duration of adiabatic half passage
%      RFe_2 - isodelay of excitation pulse
%      RFr   - duration of adiabatic full passage
%      RFr1  - start time of 1st AFP
%      RFe2  - start of 2nd AHP
%      All timings in ms
%
% Written by Joseph G. Woods, CFMRI, UCSD, March 2021
% Updated by Joseph G. Woods, University of Oxford, April 2022

function [B1, GLabel, GCont, T] = genBIR4(T, bSection, bvelCompCont)

if ~exist('T'       ,'var') || isempty(T);        error('T must be specified!'); end
if ~exist('bSection','var') || isempty(bSection); bSection = 'all'; end

% Initialise outputs in case they are not set
B1     = [];
GLabel = [];
GCont  = [];

% Set the BIR parameters (Guo and Wong, MRM 2012. http://doi.wiley.com/10.1002/mrm.24145)
wmax = 42520.0; % max frequency sweep (hz)
zeta = 43.58;   % (s^-1)
tkap = 69.65;   % tan of kappa

% Set the gradient polarities
T.polLabel = [ 1, 1];
if bvelCompCont; error('Cannot use velocity compensated control with BIR-4!');
else;            T.polCont = [ 0, 0]; end

%% Generate the BIR AHP pulses

if strcmpi(bSection,'excite') || strcmpi(bSection,'all')

    [rho, theta, ~] = genBIR(wmax, zeta, tkap, T.RFe, T.RFUP);
    
    T.B1_excite1 = T.B1max * rho.birleft  .* exp(1i * theta.birleft );
    T.B1_excite2 = T.B1max * rho.birright .* exp(1i * theta.birright);
    
    % No gradients on
    T.grad_excite1 = zeros(size(T.B1_excite1));
    T.grad_excite2 = zeros(size(T.B1_excite2));

end

%% Generate the BIR AFP pulses

if strcmpi(bSection,'refocus') || strcmpi(bSection,'all')
    
    if ~exist('rho','var')
        [rho, theta] = genBIR(wmax, zeta, tkap, T.RFe, T.RFUP);
    end
    
    T.B1_refocus = T.B1max * rho.birmid .* exp(1i * theta.birmid);

end

%% Generate the velocity encoding gradients

if strcmpi(bSection,'VSgrad') || strcmpi(bSection,'all')
    
    gLabel = genVSGrad(T, T.polLabel);
    gCont  = genVSGrad(T, T.polCont );
    
    % Increase resolution of waveform to match the RF pulses
    T.gLabel = increaseres(gLabel, round(T.GUP/T.RFUP));
    T.gCont  = increaseres(gCont , round(T.GUP/T.RFUP));
    
end
    
%% Combine the whole VS module
    
if strcmpi(bSection,'combine') || strcmpi(bSection,'all')
    
    % Gaps for VS gradients
    gap1 = zeros(round(T.RFr1*1e3/T.RFUP), 1);                % Gap during G1
    gap2 = zeros(round((T.RFe2-T.RFr1-T.RFr)*1e3/T.RFUP), 1); % Gap during G2

    %    [     excite ;  G1 ;      refocus ;  G2 ; excite       ];
    B1 = [T.B1_excite1; gap1; T.B1_refocus ; gap2; T.B1_excite2 ];
    
    %        [       excite ; VS grad ;        excite ]
    GLabel = [T.grad_excite1; T.gLabel; T.grad_excite2];
    GCont  = [T.grad_excite1; T.gCont ; T.grad_excite2];
    
end
 
end