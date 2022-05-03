%% Function to generate an Fourier-transform velocity selective inversion velocity selective module
%
% [B1, gTag, gCont, T] = genFTVSI(T, bSection)
%
% in:
%      T          - struct containing velocity selective module timings (ms
%                   except GUP and RFUP in µs).
%      Section    - which section of the module to generate: 'excite',
%                   'refocus', 'VSgrad', 'combine', or 'FTVSI'.
%                 - option 'FTVSI' generates and combines all sections.
%                 - the options can be combined to generate multiple
%                   sections at once, e.g. 'VSgradcombine' generates the
%                   VSgrad section and combines VSgrad with premade
%                   excitation and refocussing sections passed in with T.
%      bcomposite - true for composite refocussing pulses (90x-180y-90x),
%                   false for hard refocussing pulses (180x).
%
% out:
%      B1 - complex B1+ waveform (Gauss - default)
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
%      RFe2   - start time of last excitation pulse
%      RFr    - duration of refocussing pulse
%      RFr1   - start time of 1st refocussing pulse
%      RFr2   - start time of 2nd refocussing pulse
%      RFrpad - time between the HS pulses
%      All timings in ms
%
% Written by Joseph G. Woods, CFMRI, UCSD, June 2020

function [B1, gTag, gCont, T] = genFTVSI(T, section, bvelCompCont, bcomposite)

if ~exist('T'         ,'var') || isempty(T);          error('T must be specified!'); end
if ~exist('section'   ,'var') || isempty(section);    section    = 'all'; end
if ~exist('bcomposite','var') || isempty(bcomposite); bcomposite = true;    end

% Initialise outputs in case they are not set
B1    = [];
gTag  = [];
gCont = [];

if ~isfield(T,'gamrad')
    switch T.units
        case 'G' ; T.gamrad = gyroratio('rad/s/G');
        case 'T' ; T.gamrad = gyroratio('rad/s/T');
        case 'Hz'; T.gamrad = 2*pi; % Simply do not divide by γ in B1max
        otherwise; error('Units can only be G, T, or Hz!');
    end
end

% Set the gradient polarities
T.polTag = [ 1,-1, 1,-1];
if bvelCompCont; T.polCont = [ 1, 1, 1, 1];
else;            T.polCont = [ 0, 0, 0, 0]; end

%% Generate the hard excitation pulse

if contains(section,'excite') || strcmp(section,'all')

    FA     = 180/T.Nk; % excitation flip angle (degrees)
    phaseE = 0;        % phase of excitation pulse (+x)
    
    % Generate the hard pulse
    T.RFe   = 1e-3 * ceil(1e6*FA*pi/180/(T.gamrad*T.B1max)/T.GUP)*T.GUP; % (ms)
    T.RFe_2 = T.RFe/2;                                                    % approximate isodelay (for off-resonance robustness)
    T.B1_excite = genhard(FA, phaseE, T.RFe, T.B1max, T.RFUP, T.units);

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
        T.RFr = dur1 + dur2 + dur1;                                % total duration
        B1_1  = genhard(FA/2,  0, dur1, T.B1max, T.RFUP, T.units); %  90° has 0° phase
        B1_2  = genhard(FA  , 90, dur2, T.B1max, T.RFUP, T.units); % 180° has 90° phase
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
    
    gTag_VS  = genVSGrad(T, T.polTag );
    gCont_VS = genVSGrad(T, T.polCont);
    
    % Increase resolution of waveform to match the RF pulses
    T.gTag_VS  = increaseres(gTag_VS,  round(T.GUP/T.RFUP));
    T.gCont_VS = increaseres(gCont_VS, round(T.GUP/T.RFUP));
    
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
        %      [           excite;  G1 ;      refocus           ;G2+G3;      refocus              ;  G4 ];
        B1   = [B1;   T.B1_excite; gap1; T.B1_refocus(:,pc(ii)) ; gap2; T.B1_refocus(:,pc(ii+1))  ; gap3]; % 1st phase
        B1_2 = [B1_2; T.B1_excite; gap1; T.B1_refocus2(:,pc(ii)); gap2; T.B1_refocus2(:,pc(ii+1)) ; gap3]; % 2nd phase
        B1_3 = [B1_3; T.B1_excite; gap1; T.B1_refocus3(:,pc(ii)); gap2; T.B1_refocus3(:,pc(ii+1)) ; gap3]; % 3rd phase
        B1_4 = [B1_4; T.B1_excite; gap1; T.B1_refocus4(:,pc(ii)); gap2; T.B1_refocus4(:,pc(ii+1)) ; gap3]; % 4th phase
    end
    B1     = [B1  ; T.B1_excite]; % Last excitation pulse

    % Dynamic phase cycling (Liu et al. MRM 2021. https://doi.org/10.1002/mrm.28622)
    T.B1(:,1) = B1;
    T.B1(:,2) = [B1_2; T.B1_excite];
    T.B1(:,3) = [B1_3; T.B1_excite];
    T.B1(:,4) = [B1_4; T.B1_excite];
    
    gTag = [T.grad_excite; T.gTag_VS; T.grad_excite; T.gTag_VS; ...
            T.grad_excite; T.gTag_VS; T.grad_excite; T.gTag_VS; ...
            T.grad_excite; T.gTag_VS; T.grad_excite; T.gTag_VS; ...
            T.grad_excite; T.gTag_VS; T.grad_excite; T.gTag_VS; T.grad_excite ];
     
    gCont = [T.grad_excite; T.gCont_VS; T.grad_excite; T.gCont_VS; ...
             T.grad_excite; T.gCont_VS; T.grad_excite; T.gCont_VS; ...
             T.grad_excite; T.gCont_VS; T.grad_excite; T.gCont_VS; ...
             T.grad_excite; T.gCont_VS; T.grad_excite; T.gCont_VS; T.grad_excite ];
    
end
 
end
