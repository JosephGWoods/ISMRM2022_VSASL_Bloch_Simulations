%% VSASL module Bloch simulations with flow
%
% T = VSASL_ECOpt(vsType, Vcut, B1max, Gmax, SRmax, pad1, pad2, RFUP, GUP, units, bplotVS, bvelCompCont, bcomposite, bsinc)
%
% in:
%      vsType       - VS module type ('DRHS', 'DRHT', 'BIR4', 'BIR8', 'FTVSI')
%      Vcut         - velocity cutoff (cm/s)
%      B1max        - B1+ maximum amplitdue (units)
%      Gmax         - maximum gradient amplitude (units/cm)
%      SRmax        - maximum gradient slew rate (units/cm/s)
%      pad1         - pre-gradient padding (ms)
%      pad2         - post-gradient padding (ms)
%      RFUP         - RF update time (µs)
%      GUP          - gradient update time (µs)
%      units        - 'Hz', 'T', or 'G' (Hertz, Tesla, or Gauss)
%      bplotMz      - flag to plot VS module
%      bvelCompCont - flag to use velocity-compensated control
%      bcomposite   - flag to use composite refocussing pulse (FT-VSI only)
%      bsinc        - flag to use sinc modulation (FT-VSI only)
%
% out:
%      B1     - complex B1 waveform (units)
%      GLabel - label module gradient waveform (units/cm)
%      GCont  - control module gradient waveform (units/cm)
%      T      - struct of gradient and RF parameters
%
% Written by Joseph G. Woods, University of Oxford, April 2022
% Edited by Dapeng Liu, Johns Hopkins University, May 2022

function [B1, GLabel, GCont, T] = genVSASL(vsType, Vcut, B1max, Gmax, SRmax, pad1, pad2, RFUP, GUP, units, bplotVS, bvelCompCont, bcomposite, bsinc)

% Gyromagnetic ratio
switch units
    case 'G';  gam = gyroratio('Hz/G');
    case 'T';  gam = gyroratio('Hz/T');
    case 'Hz'; gam = 1;
end

% Settings for the VSASL module design
T = struct('Vcut',Vcut,'B1max',B1max,'Gmax',Gmax,'SRmax',SRmax,'pad1',pad1,'pad2',pad2, ...
    'RFUP',RFUP,'GUP',GUP,'units',units,'gam',gam,'gamrad',2*pi*gam);

func = struct('vsType',vsType,'gradAmp','scale','bsinc',bsinc);
func.RUP_GRD_ms = @(A) round(ceil(round(round(A,12)*1e3/GUP,9))*GUP*1e-3, 3);
func.Vcut2m1    = @(x) pi/(T.gamrad*2*x);

% Get/set module specicific timings (some are specified in the "gen*.m" files)
switch vsType
    case 'DRHS'
        [~,~,~,T] = genDRHS(T, 'exciterefocus', bvelCompCont); % get rf timings
    case 'DRHT'
        [~,~,~,T] = genDRHT(T, 'excite', bvelCompCont); % get rf timings
        T.RFr     = 3; % adiabatic full passage duration (ms)
    case 'BIR4'
        T.RFe   = 1.5; % adiabatic half passage duration (ms)
        T.RFe_2 = 0;   % no iso-delay needed
        T.RFr   = 3;   % adiabatic full passage duration (ms)
    case 'BIR8'
        T.RFe   = 1.5; % adiabatic half passage duration (ms)
        T.RFe_2 = 0;   % no iso-delay needed
        T.RFr   = 3;   % adiabatic full passage duration (ms)
    case 'FTVSI'
        T.Nk      = 9; % number of excitations
        [~,~,~,T] = genFTVSI(T, 'exciterefocus', bvelCompCont, bcomposite, bsinc); % get RF timings
end

% Calculate VSASL gradient flat top duration to achieve specified Vcut
T = VSTimingsEqual(func, T);

% Generate the numerical VS module
switch vsType
    case 'DRHS' ; [B1, GLabel, GCont, T] = genDRHS( T, 'all', bvelCompCont);
    case 'DRHT' ; [B1, GLabel, GCont, T] = genDRHT( T, 'all', bvelCompCont);
    case 'BIR4' ; [B1, GLabel, GCont, T] = genBIR4( T, 'all', bvelCompCont);
    case 'BIR8' ; [B1, GLabel, GCont, T] = genBIR8( T, 'all', bvelCompCont);
    case 'FTVSI'; [B1, GLabel, GCont, T] = genFTVSI(T, 'all', bvelCompCont, bcomposite, bsinc);
end
T.t = (T.RFUP : T.RFUP : T.RFUP*length(B1)) * 1e-3; % module time array for plotting

% Plot VS module
if bplotVS

    % Specify units for y-axis
    if strcmp(units,'T')
        B1plot     = B1     * 1e6;     % µT
        GLabelplot = GLabel * 1e3*1e2; % mT/m
        GContplot  = GCont  * 1e3*1e2; % mT/m
        B1units    = 'µT';
        Gunits     = 'mT/m';
    else
        B1plot     = B1;     % units
        GLabelplot = GLabel; % units/cm
        GContplot  = GCont;  % units/cm
        B1units    = units;
        Gunits     = [units '/cm'];
    end

    figure('Name','Velocity selective ASL module','Units','normalized','Position',[0,0.2,1,0.6]);
    for ii = 1:2
        subplot(2,1,ii);
        set(gca,'FontSize',12);

        % Plot B1+
        yyaxis left; hold on;
        if ii == 1
            plot(T.t, abs(B1plot), 'LineWidth', 2);
        else
            plot(T.t, real(B1plot), 'LineWidth', 2);
            plot(T.t, imag(B1plot), 'LineWidth', 2);
        end
        lim = max(abs(B1plot)) * 1.1; ylim([-lim,lim]);
        ylabel(['B_1^+ (' B1units ')'],'FontSize',16);

        % Plot gradients
        yyaxis right; hold on;
        plot(T.t, GLabelplot,  '-', 'LineWidth', 2);
        plot(T.t, GContplot , '--', 'LineWidth', 2);
        plot([0,T.t(end)], [0,0], 'k-', 'LineWidth', 1);
        lim = max(abs(GLabelplot)) * 1.1;
        ylabel(['Gradient amplitude (' Gunits ')'],'FontSize',16); ylim([-lim,lim]);

        % Labels
        box on;
        xlim([0,T.t(end)]);
        xlabel('Time (ms)','FontSize',16);
        title([vsType ' module'],'FontSize',18);
        if ii == 1; legend('|B_1^+|','G_{Label}','G_{Control}','FontSize',16,'Location','southwest');
        else;       legend('B_1^+ real','B_1^+ imag','G_{Label}','G_{Control}','FontSize',16,'Location','southwest');
        end
    end
end

end
