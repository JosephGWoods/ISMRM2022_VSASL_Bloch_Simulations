%% Function to generate a BIR4 velocity selective module
%
% [B1, gTag, gCont, T] = genbir4(T, bSection)
%
% in:
%      T        - struct containing velocity selective module timings (ms
%                 except GUP in µs) 
%      bSection - which section of the module to generate: 'excite',
%                 'refocus', 'VSgrad', 'combine', or 'BIR4'
%               - option 'BIR4' generates and combines all sections
%               - the options can be combined to generate multiple sections
%                 at once, e.g. 'VSgradcombine' generates the VSgrad
%                 section and combines VSgrad with premade RFe and RFr
%                 sections passed in with T
%
% out:
%      B1 - complex B1+ waveform
%      g  - gradient waveform (same units as T.Gmax)
%      T  - updated struct with the generated RF and gradient waveforms
%
% T parameter descriptions:
%      GUP   - gradient update time (µs)
%      RFUP  - RF update time (µs)
%      Gmax  - Max gradient amplitude (units/cm)
%      B1max - Max RF amplitude (units)
%      units - RF and gradient units ('G', 'T', 'Hz')
%      f1    - flat top time for G1
%      f2    - flat top time for G2
%      nGrad - number of VS gradients
%      r     - rise time (keep constant)
%      ta[n] - start of trapezoid attack time for nth gradient
%      td[n] - end of trapezoid decay time for nth gradient
%      RFe       - duration of adiabatic half passage for BIR-4
%      RFr       - duration of adiabatic full passage for BIR-4
%      RFr1      - start time of 1st AFP
%      RFe2      - start of 2nd AHP
%      All timings in ms
%
% Written by Joseph G. Woods, CFMRI, UCSD, March 2021

function [B1, gTag, gCont, T] = genbir4(T, bSection)

if ~exist('T'       ,'var') || isempty(T);        error('T must be specified!'); end
if ~exist('bSection','var') || isempty(bSection); bSection = 'BIR4'; end

% Initialise outputs in case they are not set
B1    = [];
gTag  = [];
gCont = [];

% Set the BIR parameters
wmax = 42520.0; % max frequency sweep (hz)
zeta = 43.58;   % (s^-1)
tkap = 69.65;   % tan of kappa

% Set the gradient polarities
if ~isfield(T,'polTag')
    T.polTag     = [ 1, 1];
    T.polEffTag  = [ 1,-1];
    T.polCont    = [ 0, 0];
    T.polEffCont = [ 0, 0];
end

%% Generate the BIR AHP pulses

if contains(bSection,'excite') || strcmp(bSection,'BIR4')

    [rho, theta, ~] = genbir(wmax, zeta, tkap, T.RFe, T.RFUP);
    
    T.B1_excite1 = T.B1max * rho.birleft  .* exp(1i * theta.birleft );
    T.B1_excite2 = T.B1max * rho.birright .* exp(1i * theta.birright);
    
    % No gradients on
    T.grad_excite1 = zeros(size(T.B1_excite1));
    T.grad_excite2 = zeros(size(T.B1_excite2));

end

%% Generate the BIR AFP pulses

if contains(bSection,'refocus') || strcmp(bSection,'BIR4')
    
    if ~exist('rho','var')
        [rho, theta] = genbir(wmax, zeta, tkap, T.RFr, T.RFUP, flip);
    end
    
    T.B1_refocus = T.B1max * rho.birmid .* exp(1i * theta.birmid); % abs(B1) and angle(B1)

end

%% Generate the velocity encoding gradients

if contains(bSection,'VSgrad') || strcmp(bSection,'BIR4')
    
    gTag_VS  = genVSGrad(T, T.polTag );
    gCont_VS = genVSGrad(T, T.polCont);
    
    % Increase resolution of waveform to match the RF pulses
    T.gTag_VS  = increaseres(gTag_VS,  round(T.GUP/T.RFUP));
    T.gCont_VS = increaseres(gCont_VS, round(T.GUP/T.RFUP));
    
end
    
%% Combine the whole VS module
    
if contains(bSection,'combine') || strcmp(bSection,'BIR4')
    
    gap1 = zeros(round(T.RFr1*1e3/T.RFUP), 1);
    gap2 = zeros(round((T.RFe2-T.RFr1-T.RFr)*1e3/T.RFUP), 1);

    %    [     excite ;  G1 ;      refocus ;  G2 ; excite       ];
    B1 = [T.B1_excite1; gap1; T.B1_refocus ; gap2; T.B1_excite2 ];
    
    %       [       excite ;   VS grad ;        excite ]
    gTag  = [T.grad_excite1; T.gTag_VS ; T.grad_excite2];
    gCont = [T.grad_excite1; T.gCont_VS; T.grad_excite2];
    
end
 
end