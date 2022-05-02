%% Simulate velocity selective modules

% Add containing folder to MATLAB path
filePath = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(filePath));

%% General settings

Vcut    = 2;             % velocity cutoff (cm/s) % CHANGE FTVSI M1 CALC SO THAT VCUT=2 CROSSES 1 AT 2 CM/S
B1max   = 20 *1e-6;      % max B1+                (µT    -> T)
Gmax    = 30 *1e-3*1e-2; % max gradient amplitude (mT/m  -> T/cm)
SRmax   = 100*1e-2;      % max gradient slew rate (T/m/s -> T/cm/s)
vspad1  = 0.1;           % time before gradients (ms)
vspad2  = 0.1;           % time after gradients (ms)
RFUP    = 10;             % RF update time (µs)
GUP     = 10;             % gradient update time (µs)
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

% Convert to Hz and Hz/cm for Bloch simulation
B1_Hz        = B1     * T.gam;
GLabel_Hzcm  = GLabel * T.gam;
GCont_Hzcm   = GCont  * T.gam;

%% Bloch simulations - velocity profile

% Simulate VSASL label and control module for a range of velocities
% Laminar flow integration is then performed

% Task: What does the velocity profile look like for different VSASL modules?
%       What happens when Vcut is changed?

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
[mzLabel_laminar, ~   ] = laminarFlowInt(mzLabel , dv);
[mzCont_laminar , ~   ] = laminarFlowInt(mzCont, dv);
[dmz_laminar    , lind] = laminarFlowInt(dmz   , dv);

% Plot results
plot_VSprofile(Vcut,dv,lind,mzLabel,mzCont,dmz,mzLabel_laminar,mzCont_laminar,dmz_laminar);

%% Bloch simulations - B1max vs velocity

% Simulates label and control module for a range of velocities and B1+
% scaling

% Some VSASL modules are more robust to B1-variation because they only use
% adiabatic pulses.

% Task: Which VSASL module(s) are most robust to B1 variation?

B1scale = linspace(0.6,1.4,101); % B1-variation (fraction)
df      = 0;                     % off-resonance (Hz)
dp      = 0;                     % position (cm)
dv      = linspace(-40,40,101);  % velocity (cm/s)

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
            'color',{'L16'},'colorbar',{'on'},'colorbarlabel',cbl,...
            'clim',{[-1,1],[-1,1],[-2,2]},'view',2,'order','columns');


%% Bloch simulations - off-resonance vs velocity

% Simulate label and control module for a range of velocities and
% off-resonances

% All of the above modules use refocussing pulses, so should be robust
% against off-resonance when B1 is as expected.

% Tasks: Does the off-resonance robustness decrease when B1scale is
%        decreased? How does this compare for different modules?

B1scale = 1;                      % B1-variation (fraction)
df      = linspace(-500,500,101); % off-resonance (Hz)
dp      = 0;                      % position (cm)
dv      = linspace(-40,40,101);   % velocity (cm/s)

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
            'color',{'L16'},'colorbar',{'on'},'colorbarlabel',cbl,...
            'clim',{[-1,1],[-1,1],[-2,2]},'view',2,'order','columns');

%% Bloch simulations - eddy currents vs Position

% Simulate label and control module for a range of positions and
% eddy-current time constants

% Eddy currents can be modelled as a convolution between the intended
% gradient and an exponential decaying gradient field, known as the eddy
% current field. The exponential decay has a time constant and amplitude.

% Eddy current fields can be non-spatially varying (B0 term) or have linear
% and even higher order components. A 'Z'-gradient can induce eddy current
% terms that vary spatially along Z or other directions, but here we just
% consider linear eddy currents terms parallel to the intended gradient.

% Tasks: Are some VSASL modules more sensitive to eddy currents?
%        What happens when vspad2 is increased?
%        What happens when a velocity-compensated control module is used?

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
            'color',{'L16'},'colorbar',{'on'},'colorbarlabel',cbl,...
            'yscale',{'log'},'clim',clim,'view',2,'order','columns');

