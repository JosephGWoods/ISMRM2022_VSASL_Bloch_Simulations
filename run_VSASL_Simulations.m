%% Simulate velocity selective modules

% Add containing folder to MATLAB path
filePath = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(filePath));

%% General settings

Vcut    = 2;             % velocity cutoff (cm/s) % CHANGE FTVSI M1 CALC SO THAT VCUT=2 CROSSES 1 AT 2 CM/S
B1max   = 20 *1e-6;      % max B1+                (µT    -> T)
Gmax    = 50 *1e-3*1e-2; % max gradient amplitude (mT/m  -> T/cm)
SRmax   = 150*1e-2;      % max gradient slew rate (T/m/s -> T/cm/s)
vspad1  = 0.1;           % padding time before gradients (ms)
vspad2  = 0.1;           % padding time after gradients (ms)
RFUP    = 10;            % RF update time (µs)
GUP     = 10;            % gradient update time (µs)
units   = 'T';           % Tesla
bplotVS = true;          % plot modules
T1      = inf;           % longitudinal relaxation time (s)
T2      = inf;           % transverse relaxation time (s)

%% Double-Refocussed Hyperbolic-Secant (DRHS) module

vsType = 'DRHS';
% vsType = 'DRHT';
% vsType = 'BIR4';
% vsType = 'BIR8';
% vsType = 'FTVSI';

% ADD options: velocity-compensated control
% ADD FTVSI options: sinc, composite refocussing

[B1, GLabel, GCont, T] = genVSASL(vsType, Vcut, B1max, Gmax, SRmax, vspad1, vspad2, RFUP, GUP, units, bplotVS);

% Convert to Hz and Hz/cm for Bloch simulations
B1_Hz        = B1     * T.gam;
GLabel_Hzcm  = GLabel * T.gam;
GCont_Hzcm   = GCont  * T.gam;

%% Bloch simulations - velocity profile

% Simulates label and control module for a range of velocities. Laminar
% flow integration is then performed.

% Task: 1. Run the simulation for each VSASL module. What does the velocity
%          profile look like for each module?
%       2. What happens when Vcut is changed?

df = 0;                     % off-resonance (Hz)
dp = 0;                     % position (cm)
dv = linspace(-40,40,1001); % velocity (cm/s)

% Run Bloch simulations
[~,~,mzLabel] = bloch_Hz(B1_Hz, GLabel_Hzcm, RFUP*1e-6, T1, T2, df, dp, dv, 0); % Label
[~,~,mzCont]  = bloch_Hz(B1_Hz, GCont_Hzcm , RFUP*1e-6, T1, T2, df, dp, dv, 0); % Control

% ASL subtraction
if contains(vsType,'FTVSI'); dmz = mzLabel - mzCont;
else;                        dmz = mzCont - mzLabel; end

% Laminar flow integration
[mzLabel_laminar, ~   ] = laminarFlowInt(mzLabel, dv);
[mzCont_laminar , ~   ] = laminarFlowInt(mzCont , dv);
[dmz_laminar    , lind] = laminarFlowInt(dmz    , dv);

% Plot results
plot_VSprofile(Vcut,dv,lind,mzLabel,mzCont,dmz,mzLabel_laminar,mzCont_laminar,dmz_laminar);

%% Bloch simulations - B1 vs velocity

% Simulates label and control module for a range of velocities and B1
% scaling

% Some VSASL modules are more robust to B1-variation because they use
% adiabatic pulses.
%  - DRHS/DRHT modules use hard excitation pulses but adiabatic refocussing
%    pulses.
%  - BIR-4/BIR-8 use both adiabatic excitation and refocussing pulses.
%  - FT-VSI uses hard excitation pulses but composite refocussing pulses

% Tasks: 1. Run the simulation for each VSASL module. Which are least/most
%           robust to B1 variation?
%        2. Are both the label or control module sensitive to B1?
%        3. For FT-VSI, try switching off the composite refocussing pulses
%           (which uses hard refocussing pulses instead). What happens?

B1scale = linspace(0.6,1.4,101); % B1-variation (fraction)
df      = 0;                     % off-resonance (Hz)
dp      = 0;                     % position (cm)
dv      = linspace(-20,20,101);  % velocity (cm/s)

% Run Bloch simulations
mzLabel = zeros(length(B1scale), length(dv));
mzCont  = zeros(length(B1scale), length(dv));
for ii = 1:length(B1scale)
    [~,~,mzLabel(ii,:)] = bloch_Hz(B1scale(ii)*B1_Hz, GLabel_Hzcm, RFUP*1e-6, T1, T2, df, dp, dv, 0); % Label
    [~,~,mzCont(ii,:)]  = bloch_Hz(B1scale(ii)*B1_Hz, GCont_Hzcm , RFUP*1e-6, T1, T2, df, dp, dv, 0); % Control
end

% ASL subtraction
if contains(vsType,'FTVSI'); dmz = mzLabel  - mzCont;
else;                        dmz = mzCont - mzLabel; end

% Plot data
data1 = {dv,dv,dv};
data2 = {B1scale,B1scale,B1scale};
data3 = {mzCont,mzLabel,dmz};
ti    = {'Control','Label','Difference'};
cbl   = {'Mz/M_0','Mz/M_0','ΔMz/M_0'};

surf_custom('data1',data1,'data2',data2,'data3',data3,...
            'xlabel',{'Velocity (cm/s)'},'ylabel',{'B1^+ scale'},'title',ti,...
            'color',{'jet'},'colorbarlabel',cbl,...
            'clim',{[-1,1],[-1,1],[0,2]});


%% Bloch simulations - off-resonance (B0) vs velocity

% Simulates label and control module for a range of velocities and
% off-resonances

% All of the above modules use refocussing pulses, so should be fairly
% robust to off-resonance variation when B1scale = 1 (nominal B1).

% Tasks: 1. Run the simulation for each VSASL module to confirm they are
%           all robust to off-resonance.
%        2. Does the off-resonance robustness decrease when B1scale is
%           decreased? How does this compare for different modules?

B1scale = 1;                      % B1-variation (fraction)
df      = linspace(-500,500,101); % off-resonance (Hz)
dp      = 0;                      % position (cm)
dv      = linspace(-20,20,101);   % velocity (cm/s)

% Run Bloch simulations
[~,~,mzLabel] = bloch_Hz(B1scale*B1_Hz, GLabel_Hzcm, RFUP*1e-6, T1, T2, df, dp, dv, 0); % Label
[~,~,mzCont ] = bloch_Hz(B1scale*B1_Hz, GCont_Hzcm , RFUP*1e-6, T1, T2, df, dp, dv, 0); % Control

% ASL subtraction
if contains(vsType,'FTVSI'); dmz = mzLabel - mzCont;
else;                        dmz = mzCont - mzLabel; end

% Plot data
data1 = {dv,dv,dv};
data2 = {df,df,df};
data3 = {mzCont,mzLabel,dmz};
ti    = {'Control','Label','Difference'};
cbl   = {'Mz/M_0','Mz/M_0','ΔMz/M_0'};

surf_custom('data1',data1,'data2',data2,'data3',data3,...
            'xlabel',{'Velocity (cm/s)'},'ylabel',{'Off-resonance (Hz)'},'title',ti,...
            'color',{'jet'},'colorbarlabel',cbl,...
            'clim',{[-1,1],[-1,1],[0,2]});

%% Bloch simulations - eddy currents vs Position

% Simulates label and control module for a range of positions and
% eddy-current time constants

% Eddy currents can be modelled as a convolution between the intended
% gradient and an exponential decaying gradient field, known as the eddy
% current field. The exponential decay has a time constant and amplitude.

% Eddy current fields can be non-spatially varying (B0 term) or have linear
% or higher order components. A Z-gradient can induce eddy current fields
% that vary spatially along Z or other directions, but here we just
% consider linear eddy currents terms parallel to the intended gradient.

% Tasks: 1. Run the simulation for each VSASL module. Are some VSASL
%           modules more sensitive to eddy currents than others? 
%        2. What happens when vspad2 is increased to 2 ms?
%        3. What happens when a velocity-compensated control module is used?
%        4. (ADVANCED) Alter code to plot eddy currents against velocity
%           for a fixed position of 24 cm to demonstrate how the velocity
%           profile is distorted by eddy current effects.

% Generate linear eddy currents parallel to VS gradients
tau = logspace(3,-1,100); % eddy current time constants (ms)
A   = 0.0025;             % eddy current amplitude (fraction of input gradient)
GLabel_EC_Hzcm = GLabel_Hzcm + gradec(GLabel_Hzcm, A, tau, GUP); % add eddy currents to gradients
GCont_EC_Hzcm  = GCont_Hzcm  + gradec(GCont_Hzcm , A, tau, GUP); % add eddy currents to gradients

B1scale = 1;                    % B1-variation (fraction)
df      = 0;                    % off-resonance (Hz)
dp      = linspace(-24,24,101); % position (cm)
dv      = 0;                    % velocity (cm/s)

% Run Bloch simulations
mzLabel = zeros(length(tau),length(dp));
mzCont  = zeros(length(tau),length(dp));
for ii = 1:length(tau)
    [~,~,mzLabel(ii,:)] = bloch_Hz(B1scale*B1_Hz, GLabel_EC_Hzcm(:,ii), RFUP*1e-6, T1, T2, df, dp, dv, 0); % Label
    [~,~,mzCont(ii,:)]  = bloch_Hz(B1scale*B1_Hz, GCont_EC_Hzcm(:,ii) , RFUP*1e-6, T1, T2, df, dp, dv, 0); % Control
end

% ASL subtraction
if contains(vsType,'FTVSI'); dmz = mzLabel - mzCont;
else;                        dmz = mzCont - mzLabel; end

% Plot data
data1 = {dp,dp,dp};
data2 = {tau,tau,tau};
data3 = {mzCont,mzLabel,dmz};
ti    = {'Control','Label','Difference'};
cbl   = {'Mz/M_0','Mz/M_0','ΔMz/M_0'};
if contains(vsType,'FTVSI'); clim = {[-1,0],[-1,0],[0,1]};
else;                        clim = {[ 0,1],[ 0,1],[0,1]}; end

surf_custom('data1',data1,'data2',data2,'data3',data3,...
            'xlabel',{'Position (cm)'},'ylabel',{'Time constant (ms)'},'title',ti,...
            'color',{'jet'},'colorbarlabel',cbl,...
            'yscale',{'log'},'clim',clim);

%% Bloch simulations - B1 vs position

% Simulates label and control module for a range of positions and
% B1 scaling

% Devations from the nominal B1 can lead to inaccurate refocussing for
% FT-VSI (much less so for other VSASL modules that use adiabatic
% refocussing). This means that the phase accumulation by the
% velocity-encoding gradients is not perfectly reversed for static spins,
% creating spatial "stripe"-artifacts. (see Liu et al. MRM 2021.
% https://doi.org/10.1002/mrm.28622)
% 
% Liu et al. proposed using dynamic phase cycling, where a 90° phase is
% added to the refocussing pulses over x4 TRs. This spatially shifts the
% stripe artifact in each acquisition so that it cancels out after
% averaging.
%
% This simulation demonstrates this effect and solution.

% Tasks: 1. Run the simulation for using FT-VSI.
%        2. Does phase cycling still work when off-resonance = 200 Hz?
%        3. Is this a strong effect for the other VSASL modules?

B1scale = linspace(0.6,1.4,101);  % B1-variation (fraction)
df      = 0;                      % off-resonance (Hz)
dp      = linspace(-0.2,0.2,101); % position (cm)
dv      = 0;                      % velocity (cm/s)

% Check not BIR-4
if strcmpi(vsType,'BIR4')
    error('Dynamic phase cycling not compatible with BIR-4!');
end

dynphase = size(T.B1,2); % x4 TR phase cycling
mzLabel  = cell(1,dynphase);
mzCont   = cell(1,dynphase);
dmz      = cell(1,dynphase);
for jj = 1:dynphase

    % Update B1 waveform (convert to Hz)
    B1_Hz = T.B1(:,jj) * T.gam; % phase cycle refocussing pulses

    % Run Bloch simulations
    for ii = 1:length(B1scale)
        [~,~,mzLabel{jj}(ii,:)] = bloch_Hz(B1scale(ii)*B1_Hz, GLabel_Hzcm, RFUP*1e-6, T1, T2, df, dp, dv, 0);
        [~,~,mzCont{jj}(ii,:)]  = bloch_Hz(B1scale(ii)*B1_Hz, GCont_Hzcm , RFUP*1e-6, T1, T2, df, dp, dv, 0);
    end

    % ASL subtraction
    if contains(vsType,'FTVSI'); dmz{jj} = mzLabel{jj} - mzCont{jj};
    else;                        dmz{jj} = mzCont{jj} - mzLabel{jj}; end

end

% Mean over phase cycled TRs
mzLabelMean = mean(reshape(cell2mat(mzLabel),length(B1scale),length(dp),dynphase),3);
mzContMean  = mean(reshape(cell2mat(mzCont) ,length(B1scale),length(dp),dynphase),3);
dmzMean     = mean(reshape(cell2mat(dmz)    ,length(B1scale),length(dp),dynphase),3);

% Plot data
data1 = repmat({dp,dp,dp},1,dynphase+1);
data2 = repmat({B1scale,B1scale,B1scale},1,dynphase+1);
data3 = [mzCont,mzContMean,mzLabel,mzLabelMean,dmz,dmzMean];
ti    = [repmat({'Control'},1,dynphase),{'Control mean'},...
         repmat({'Label'}  ,1,dynphase),{'LabeL mean'}   ,...
         repmat({'Difference'},1,dynphase),{'Difference mean'}];
cbl   = [repmat({'Mz/M_0'},1,2*(dynphase+1)),repmat({'ΔMz/M_0'},1,dynphase+1)];
clim  = [repmat({[-1,0]},1,2*(dynphase+1)),repmat({[-0.5,0.5]},1,dynphase+1)];

surf_custom('data1',data1,'data2',data2,'data3',data3,...
            'xlabel',{'Position (cm)'},'ylabel',{'B1^+ scale'},'title',ti,...
            'color',{'jet'},'colorbarlabel',cbl,...
            'clim',clim,'dim',[3,dynphase+1]);
