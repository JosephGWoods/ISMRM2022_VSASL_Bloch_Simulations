%% Function to generate a double-refocussed hyperbolic-tangent velocity selective module
%
% [B1, gTag, gCont, T] = genDRHT(T, bSection)
%
% in:
%      T        - struct containing velocity selective module timings (ms
%                 except GUP in µs) 
%      bSection - which section of the module to generate: 'excite',
%                 'refocus', 'VSgrad', 'combine', or 'DRHT'
%               - option 'DRHT' generates and combines all sections
%               - the options can be combined to generate multiple
%                 sections at once, e.g. 'VSgradcombine' generates the
%                 VSgrad section and combines VSgrad with premade
%                 excitation and refocussing sections passed in with T.
%
% out:
%      B1 - complex B1+ waveform
%      g  - gradient waveform (same units as T.Gmax)
%      T  - updated struct with the generated RF and gradient waveforms
%
% T parameter descriptions:
%      GUP    - gradient update time (µs)
%      RFUP   - RF update time (µs)
%      Gmax   - max gradient amplitude (units/cm)
%      SRmax  - max gradient slew rate (units/cm/s
%      B1max  - max B1+ amplitude (units)
%      units  - B1+ and gradient units ('G', 'T', 'Hz')
%      f      - gradient flat top time
%      r      - gradient rise time
%      ta[n]  - start of trapezoid attack time for nth gradient
%      td[n]  - end of trapezoid decay time for nth gradient
%      RFe    - duration of excitation pulse
%      RFe_2  - isodelay of excitation pulse
%      RF2    - start time of last excitation pulse
%      RFr    - duration of refocussing pulse
%      RFr1   - start time of 1st refocussing pulse
%      RFr2   - start time of 2nd refocussing pulse
%      RFrpad - time between the HS pulses
%      All timings in ms
%
% Written by Joseph G. Woods, University of Oxford, April 2022

function [B1, gTag, gCont, T] = genDRHT(T, bSection)

if ~exist('T'       ,'var') || isempty(T);        error('T must be specified!'); end
if ~exist('bSection','var') || isempty(bSection); bSection = 'DRHT'; end

% Initialise outputs in case they are not set
B1    = [];
gTag  = [];
gCont = [];

RUP_GRD_ms = @(A) round(ceil(round(round(A,12)*1e3/T.GUP,9))*T.GUP*1e-3, 3);

% Set the gradient polarities
if ~isfield(T,'polTag')
    T.polTag     = [ 1,-1, 1,-1];
    T.polEffTag  = [ 1, 1,-1,-1];
    T.polCont    = [ 0, 0, 0, 0];
    T.polEffCont = [ 0, 0, 0, 0];
    %T.PolCont    = [ 1, 1, 1, 1];
    %T.PolEffCont = [ 1,-1,-1, 1];
    %T.PolCont    = [-1,-1,-1,-1];
    %T.PolEffCont = [-1, 1, 1,-1];
end

%% Generate the hard excitation pulses

if contains(bSection,'excite') || strcmp(bSection,'DRHT')

    FA     = 90;  % Flip angle (degrees)
    phase1 = 0;   % Phase of flip down
    phase2 = 180; % Phase of flip up
    
    % Generate the hard pulses
    T.RFe = 1e-3 * ceil(1e6*FA*pi/180/(T.gamrad*T.B1max)/T.GUP)*T.GUP; % ms
    T.B1_excite1 = genhard(FA, phase1, [], T.B1max, T.RFUP, T.units);
    T.B1_excite2 = genhard(FA, phase2, [], T.B1max, T.RFUP, T.units);
    
    % Approximate isodelay (for off-resonance robustness)
    vsisd = 1.3125;
    T.RFe_2 = T.RFe - RUP_GRD_ms(T.RFe/2/vsisd);

    % No gradients on
    T.grad_excite1 = zeros(size(T.B1_excite1));
    T.grad_excite2 = zeros(size(T.B1_excite2));

end

%% Generate the HT refocussing pulses
% Uses hyperbolic-tanget pulses (see Qin and van Zijl. MRM 2016 https://doi.org/10.1002/mrm.26010)

if contains(bSection,'refocus') || strcmp(bSection,'DRHT')

    % Set the BIR parameters (Guo and Wong, MRM 2012. http://doi.wiley.com/10.1002/mrm.24145)
    wmax = 42520.0;   % max frequency sweep (hz)
    zeta = 43.58;     % (s^-1)
    tkap = 69.65;     % tan of kapp
    phase = [-90,90]; % Initial phase of refocussing pulses (180° phase between them)
    
    [rho, theta, ~] = genbir(wmax, zeta, tkap, RUP_GRD_ms(T.RFr/2), T.RFUP);

    T.B1_refocus(:,1) = T.B1max * rho.birmid .* exp(1i * theta.birmid) * exp(1i * phase(1)*pi/180); % abs(B1) and angle(B1)
    T.B1_refocus(:,2) = T.B1max * rho.birmid .* exp(1i * theta.birmid) * exp(1i * phase(2)*pi/180); % abs(B1) and angle(B1)
    
    % Generate the 90° phase increments (Liu et al. MRM 2021. https://doi.org/10.1002/mrm.28622)
    T.B1_refocus2 = T.B1_refocus * exp(1i *  90*pi/180);
    T.B1_refocus3 = T.B1_refocus * exp(1i * 180*pi/180);
    T.B1_refocus4 = T.B1_refocus * exp(1i * 270*pi/180);
    
end

%% Generate the velocity encoding gradients

if contains(bSection,'VSgrad') || strcmp(bSection,'DRHT')
    
    gTag_VS  = genVSGrad(T, T.polTag );
    gCont_VS = genVSGrad(T, T.polCont);
    
    % Increase resolution of waveform to match the RF pulses
    T.gTag_VS  = increaseres(gTag_VS,  round(T.GUP/T.RFUP));
    T.gCont_VS = increaseres(gCont_VS, round(T.GUP/T.RFUP));
    
end
    
%% Combine the whole VS module
    
if contains(bSection,'combine') || strcmp(bSection,'DRHT')

    gap1 = zeros(round(T.RFr1*1e3/T.RFUP), 1);
    gap2 = zeros(round((T.RFr2-T.RFr1-T.RFr)*1e3/T.RFUP), 1);
    gap3 = zeros(round((T.RFe2-T.RFr2-T.RFr)*1e3/T.RFUP), 1);
    
    %        [     excite;   G1 ;      refocus ;G2+G3;      refocus ;  G4 ;      excite ];
    B1     = [T.B1_excite1; gap1; T.B1_refocus(:,1) ; gap2; T.B1_refocus(:,2) ; gap3; T.B1_excite2];
    T.B1_2 = [T.B1_excite1; gap1; T.B1_refocus2(:,1); gap2; T.B1_refocus2(:,2); gap3; T.B1_excite2]; % 2nd phase
    T.B1_3 = [T.B1_excite1; gap1; T.B1_refocus3(:,1); gap2; T.B1_refocus3(:,2); gap3; T.B1_excite2]; % 3rd phase 
    T.B1_4 = [T.B1_excite1; gap1; T.B1_refocus4(:,1); gap2; T.B1_refocus4(:,2); gap3; T.B1_excite2]; % 4th phase
    
    gTag  = [T.grad_excite1; T.gTag_VS ; T.grad_excite2];
    gCont = [T.grad_excite1; T.gCont_VS; T.grad_excite2];
    
end
 
end
