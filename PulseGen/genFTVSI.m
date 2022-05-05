%% Function to generate an Fourier-transform velocity selective inversion velocity selective module
%
% [B1, GLabel, GCont, T] = genFTVSI(T, bSection)
%
% in:
%      T            - struct of gradient and RF parameters 
%      bSection     - sections of the module to generate: 'excite',
%                     'refocus', 'VSgrad', 'combine', or 'all'
%      bvelCompCont - flag to use velocity-compensated control
%      bcomposite   - flag to use composite refocussing pulses
%                     (90x-180y-90x). If false, hard refocussing pulses
%                     (180x) are used
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
%      f      - gradient flat top time
%      r      - gradient rise time
%      ta[n]  - start of trapezoid attack time for nth gradient
%      td[n]  - end of trapezoid decay time for nth gradient
%      RFe    - duration of excitation pulse
%      RFe_2  - isodelay of excitation pulse
%      RFe2   - start time of last excitation pulse
%      RFr    - duration of refocussing pulse
%      RFr1   - start time of 1st refocussing pulse
%      RFr2   - start time of 2nd refocussing pulse
%      RFrpad - time between the HS pulses
%      Nk     - number of FT-VSI excitation pulses
%      All timings in ms
%
% Written by Joseph G. Woods, CFMRI, UCSD, June 2020
% Updated by Joseph G. Woods, University of Oxford, April 2022
% Edited by Dapeng Liu, Johns Hopkins University, May 2022

function [B1, GLabel, GCont, T] = genFTVSI(T, section, bvelCompCont, bcomposite, bsinc)

if ~exist('T'         ,'var') || isempty(T);          error('T must be specified!'); end
if ~exist('section'   ,'var') || isempty(section);    section    = 'all'; end
if ~exist('bcomposite','var') || isempty(bcomposite); bcomposite = true;    end

% Initialise outputs in case they are not set
B1     = [];
GLabel = [];
GCont  = [];

if ~isfield(T,'gamrad')
    switch T.units
        case 'G' ; T.gamrad = gyroratio('rad/s/G');
        case 'T' ; T.gamrad = gyroratio('rad/s/T');
        case 'Hz'; T.gamrad = 2*pi; % Simply do not divide by γ in B1max
        otherwise; error('Units can only be G, T, or Hz!');
    end
end

% Set the gradient polarities
T.polLabel = [ 1,-1, 1,-1];
if bvelCompCont; T.polCont = [ 1, 1, 1, 1];
else;            T.polCont = [ 0, 0, 0, 0]; end

%% Generate the hard excitation pulse

if contains(section,'excite') || strcmp(section,'all')

    FA     = 180; % total excitation flip angle (degrees)
    phaseE = 0;   % phase of excitation pulse (+x)

    if bsinc
        % Sinc modulation (Guo et al. MRM 2021 https://doi.org/10.1002/mrm.28572)
        FTmod = sinc((1:T.Nk)/ceil(T.Nk/2)-1);                                        % generate Nk-point single-lobe non-windowed sinc
        FTmod = T.Nk * FTmod / sum(FTmod);                                            % normalise sinc modulation by T.Nk
        FA    = FA/T.Nk;                                                              % mean excitation flip angle (degrees)
        T.RFe = 1e-3 * ceil(1e6*max(FTmod)*FA*pi/180/(T.gamrad*T.B1max)/T.GUP)*T.GUP; % excitation duration (ms)
    else
        FA    = FA/T.Nk;                                                   % excitation flip angle (degrees)
        T.RFe = 1e-3 * ceil(1e6*FA*pi/180/(T.gamrad*T.B1max)/T.GUP)*T.GUP; % excitation duration (ms)
        FTmod = ones(T.Nk,1);                                              % rectangular modulation (i.e. no modulation)
    end
    
    % Generate the hard pulse
    T.RFe_2     = T.RFe/2;                                              % approximate isodelay (for off-resonance robustness)
    T.B1_excite = genhard(FA, phaseE, T.RFe, T.B1max, T.RFUP, T.units); % 

    % No gradients on
    T.grad_excite = zeros(size(T.B1_excite));

end

%% Generate the hard refocussing pulse

if contains(section,'refocus') || strcmp(section,'all')

    FA     = 180;      % Flip angle (degrees)
    phaseR = [90,-90]; % Phase of refocussing pulses [+y, -y]
    
    if ~bcomposite
        % Generate the hard pulse
        T.RFr = 1e-3 * ceil(1e6*FA*pi/180/(T.gamrad*T.B1max)/T.GUP)*T.GUP; % (ms)
        B1_refocus = genhard(FA, 0, T.RFr, T.B1max, T.RFUP, T.units);
    else
        % Generate the composite pulse (see Liu et al. MRM 2021. doi:https://doi.org/10.1002/mrm.28310)
        dur1  = 1e-3 * ceil(1e6*(FA/2)*pi/180/(T.gamrad*T.B1max)/T.GUP)*T.GUP; %  90° (ms)
        dur2  = 1e-3 * ceil(1e6* FA   *pi/180/(T.gamrad*T.B1max)/T.GUP)*T.GUP; % 180° (ms)
        T.RFr = dur1 + dur2 + dur1;                                            % total duration
        B1_1  = genhard(FA/2,  0, dur1, T.B1max, T.RFUP, T.units);             %  90° has 0° phase
        B1_2  = genhard(FA  , 90, dur2, T.B1max, T.RFUP, T.units);             % 180° has 90° phase
        B1_refocus = [B1_1; B1_2; B1_1];
    end
    
    % Add phase for +y and -y pulses (for MLEV-16)
    T.B1_refocus      = zeros(length(B1_refocus),2);
    T.B1_refocus(:,1) = B1_refocus .* exp(1i * phaseR(1)*pi/180); % +y
    T.B1_refocus(:,2) = B1_refocus .* exp(1i * phaseR(2)*pi/180); % -y
    
    % Generate the 90° phase increments (Liu et al. MRM 2021. https://doi.org/10.1002/mrm.28622)
    T.B1_refocus2 = T.B1_refocus * exp(1i *  90*pi/180);
    T.B1_refocus3 = T.B1_refocus * exp(1i * 180*pi/180);
    T.B1_refocus4 = T.B1_refocus * exp(1i * 270*pi/180);
    
    % No gradients on
    T.grad_refocus = zeros(length(T.B1_refocus),1);
    
end

%% Generate the velocity encoding gradients

if contains(section,'VSgrad') || strcmp(section,'all')
    
    gLabel = genVSGrad(T, T.polLabel);
    gCont  = genVSGrad(T, T.polCont );
    
    % Increase resolution of waveform to match the RF pulses
    T.gLabel = increaseres(gLabel, round(T.GUP/T.RFUP));
    T.gCont  = increaseres(gCont , round(T.GUP/T.RFUP));
    
end
    
%% Combine the whole VS module
    
if contains(section,'combine') || strcmp(section,'all')

    % Gaps for VS gradients
    gap1 = zeros(round( T.RFr1              *1e3/T.RFUP), 1); % Gap during G1
    gap2 = zeros(round((T.RFr2-T.RFr1-T.RFr)*1e3/T.RFUP), 1); % Gap during G2 and G3
    gap3 = zeros(round((T.RFe2-T.RFr2-T.RFr)*1e3/T.RFUP), 1); % Gap during G4
    
    pc = [1,1,2,2,2,1,1,2,2,2,1,1,1,2,2,1]; % MLEV-16 order 1

    B1 = []; B1_2 = []; B1_3 = []; B1_4 = [];
    for ii = 1 : 2 : length(pc)        
        %      [           excite                ;  G1 ;      refocus           ;G2+G3;      refocus              ;  G4 ];
        B1   = [B1;   T.B1_excite*FTmod((ii+1)/2); gap1; T.B1_refocus(:,pc(ii)) ; gap2; T.B1_refocus(:,pc(ii+1))  ; gap3]; % 1st phase
        B1_2 = [B1_2; T.B1_excite*FTmod((ii+1)/2); gap1; T.B1_refocus2(:,pc(ii)); gap2; T.B1_refocus2(:,pc(ii+1)) ; gap3]; % 2nd phase
        B1_3 = [B1_3; T.B1_excite*FTmod((ii+1)/2); gap1; T.B1_refocus3(:,pc(ii)); gap2; T.B1_refocus3(:,pc(ii+1)) ; gap3]; % 3rd phase
        B1_4 = [B1_4; T.B1_excite*FTmod((ii+1)/2); gap1; T.B1_refocus4(:,pc(ii)); gap2; T.B1_refocus4(:,pc(ii+1)) ; gap3]; % 4th phase
    end
    B1 = [B1; T.B1_excite*FTmod(end)]; % Last excitation pulse

    % Dynamic phase cycling (Liu et al. MRM 2021. https://doi.org/10.1002/mrm.28622)
    T.B1(:,1) = B1;
    T.B1(:,2) = [B1_2; T.B1_excite*FTmod(end)];
    T.B1(:,3) = [B1_3; T.B1_excite*FTmod(end)];
    T.B1(:,4) = [B1_4; T.B1_excite*FTmod(end)];
    
    GLabel = [T.grad_excite; T.gLabel; T.grad_excite; T.gLabel; ...
              T.grad_excite; T.gLabel; T.grad_excite; T.gLabel; ...
              T.grad_excite; T.gLabel; T.grad_excite; T.gLabel; ...
              T.grad_excite; T.gLabel; T.grad_excite; T.gLabel; T.grad_excite];
     
    GCont = [T.grad_excite; T.gCont; T.grad_excite; T.gCont; ...
             T.grad_excite; T.gCont; T.grad_excite; T.gCont; ...
             T.grad_excite; T.gCont; T.grad_excite; T.gCont; ...
             T.grad_excite; T.gCont; T.grad_excite; T.gCont; T.grad_excite];
    
end
 
end
