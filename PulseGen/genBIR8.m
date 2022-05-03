%% Function to generate a BIR-8 velocity selective module
%
% [B1, gTag, gCont, T] = genBIR8(T, bSection)
%
% in:
%      T        - struct containing velocity selective module timings (ms
%                 except GUP in µs) 
%      bSection - which section of the module to generate: 'excite',
%                 'refocus', 'VSgrad', 'combine', or 'BIR8'
%               - option 'BIR8' generates and combines all sections
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
%      Gmax  - max gradient amplitude (units/cm)
%      SRmax - max gradient slew rate (units/cm/s
%      B1max - max B1+ amplitude (units)
%      units - B1+ and gradient units ('G', 'T', 'Hz')
%      f     - gradient flat top time
%      r     - gradient rise time
%      ta[n] - start of trapezoid attack time for nth gradient
%      td[n] - end of trapezoid decay time for nth gradient
%      RFe   - duration of adiabatic half passage
%      RFe_2 - isodelay of excitation pulse
%      RFr   - duration of adiabatic full passage
%      RFr1  - start time of 1st AFP
%      RFr2  - start time of 2nd AFP
%      RFr3  - start time of 3rd AFP
%      RFe2  - start of 2nd AHP
%      All timings in ms
%
% Written by Joseph G. Woods, CFMRI, UCSD, June 2020

function [B1, gTag, gCont, T] = genBIR8(T, bSection)

if ~exist('T'       ,'var') || isempty(T);        error('T must be specified!'); end
if ~exist('bSection','var') || isempty(bSection); bSection = 'BIR8'; end

% Initialise outputs in case they are not set
B1    = [];
gTag  = [];
gCont = [];

% Set the BIR parameters (Guo and Wong, MRM 2012. http://doi.wiley.com/10.1002/mrm.24145)
wmax = 42520.0; % max frequency sweep (hz)
zeta = 43.58;   % (s^-1)
tkap = 69.65;   % tan of kapp

% Set the gradient polarities
if ~isfield(T,'polTag')
    T.polTag     = [ 1,-1,-1, 1];
    T.polEffTag  = [ 1, 1,-1,-1];
    T.polCont    = [ 0, 0, 0, 0];
    T.polEffCont = [ 0, 0, 0, 0];
end

%% Generate the BIR AHP pulses

if contains(bSection,'excite') || strcmp(bSection,'BIR8')

    [rho, theta, ~] = genbir(wmax, zeta, tkap, T.RFe, T.RFUP);
    
    T.B1_excite1 = T.B1max * rho.birleft  .* exp(1i * theta.birleft );
    T.B1_excite2 = T.B1max * rho.birright .* exp(1i * theta.birright);
    
    % No gradients on
    T.grad_excite1 = zeros(size(T.B1_excite1));
    T.grad_excite2 = zeros(size(T.B1_excite2));

end

%% Generate the BIR AFP pulses

if contains(bSection,'refocus') || strcmp(bSection,'BIR8')
    
    if ~exist('rho','var')
        [rho, theta] = genbir(wmax, zeta, tkap, T.RFe, T.RFUP);
    end
    
    T.B1_refocus = T.B1max * rho.birmid .* exp(1i * theta.birmid); % abs(B1) and angle(B1)
    
    % Generate the 90° phase increments (Liu et al. MRM 2021. https://doi.org/10.1002/mrm.28622)
    T.B1_refocus2 = T.B1_refocus * exp(1i *  90*pi/180);
    T.B1_refocus3 = T.B1_refocus * exp(1i * 180*pi/180);
    T.B1_refocus4 = T.B1_refocus * exp(1i * 270*pi/180);

end

%% Generate the velocity encoding gradients

if contains(bSection,'VSgrad') || strcmp(bSection,'BIR8')
    
    gTag_VS  = genVSGrad(T, T.polTag );
    gCont_VS = genVSGrad(T, T.polCont);
    
    % Increase resolution of waveform to match the RF pulses
    T.gTag_VS  = increaseres(gTag_VS,  round(T.GUP/T.RFUP));
    T.gCont_VS = increaseres(gCont_VS, round(T.GUP/T.RFUP));
    
end
    
%% Combine the whole VS module
    
if contains(bSection,'combine') || strcmp(bSection,'BIR8')
    
    % Gaps for VS gradients
    gap1 = zeros(round(T.RFr1*1e3/T.RFUP), 1);
    gap2 = zeros(round((T.RFr2-T.RFr1-T.RFr)*1e3/T.RFUP), 1);
    gap3 = zeros(round((T.RFr3-T.RFr2-T.RFr)*1e3/T.RFUP), 1);
    gap4 = zeros(round((T.RFe2-T.RFr3-T.RFr)*1e3/T.RFUP), 1);

    %        [     excite ;  G1 ;      refocus ;  G2 ;      refocus;  G3 ;      refocus ;  G4 ;      excite ];
    B1     = [T.B1_excite1; gap1; T.B1_refocus ; gap2; T.B1_refocus; gap3; T.B1_refocus ; gap4; T.B1_excite2]; % 1st phase

    % Dynamic phase cycling (Liu et al. MRM 2021. https://doi.org/10.1002/mrm.28622)
    T.B1(:,1) = B1;
    T.B1(:,2) = [T.B1_excite1; gap1; T.B1_refocus2; gap2; T.B1_refocus; gap3; T.B1_refocus2; gap4; T.B1_excite2]; % 2nd phase
    T.B1(:,3) = [T.B1_excite1; gap1; T.B1_refocus3; gap2; T.B1_refocus; gap3; T.B1_refocus3; gap4; T.B1_excite2]; % 3rd phase
    T.B1(:,4) = [T.B1_excite1; gap1; T.B1_refocus4; gap2; T.B1_refocus; gap3; T.B1_refocus4; gap4; T.B1_excite2]; % 4th phase
    
    %       [       excite ;   VS grad ;        excite ]
    gTag  = [T.grad_excite1; T.gTag_VS ; T.grad_excite2];
    gCont = [T.grad_excite1; T.gCont_VS; T.grad_excite2];
    
end
 
end