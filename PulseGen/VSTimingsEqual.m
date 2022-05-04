%% Functions to calculate gradient flat top times and generate velocity selective module timings
%
% T = VSTimingsEqual(func, T)
%
% in:
%      func - struct containing RUP_GRD_ms, vsType, and gradAmp options
%      T    - struct containing velocity selective module timings (ms
%             except GUP in µs) 
%
% out:
%      T  - updated struct with the generated RF and gradient timings
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
%      RFr3   - start time of 3rd refocussing pulse (BIR-8)
%      RFrpad - time between the HS/HT pulses (DRHS/DRHT/FTVSI)
%      Nk     - number of FT-VSI excitation pulses
%      All timings in ms
%
% Example timing layout for DRHS, DRHT, and FT-VSI:
%
%              f                                                  f
%             ___                                                ___
%            /   \                                              /   \
%          r/     \r                                          r/     \r
% RFe______/       \______RFr______         ______RFrpad______/       \______RFr______         ______RFe
%    vspad1         vspad2   vspad1\       /vspad2      vspad1         vspad2   vspad1\       /
%                                  r\     /r                                          r\     /r
%                                    \___/                                              \___/
%                                      f                                                  f
%
% Written by Jia Guo and Joseph G. Woods, CFMRI, UCSD

function T = VSTimingsEqual(func, T)
 
m1Target   = 1e12 * pi/(T.gamrad*2*T.Vcut); % 1st gradient moment to achieve Vcut
RUP_GRD_ms = func.RUP_GRD_ms;               % extract often used timing rounding function

dt  = T.GUP;                               % gradient update time
T.r = RUP_GRD_ms(abs(T.Gmax/T.SRmax*1e3)); % gradient rise/ramp time
R   = T.r*1e3/dt;                          % convert to rise time to number of gradient raster points
radicand = 0;

% This section essentially tries to solve the 0th and 1st moment gradient
% equations for the gradient flat top time which is the only unknown. Gmax
% is then adjusted so that Vcut is exactly achieved.
%
% i.e. solve the following equations for f (flat top time):
% m0 = 0
% m1 = π / (2*π*γ*2*Vcut)

% TODO: add comments and equations for clarity.

if strcmpi(func.vsType,'DRHS') || strcmpi(func.vsType,'DRHT')
    p        = (4*(T.vspad1+T.vspad2)+2*T.RFr+4*T.RFe_2+8*T.r)*1e3/dt; % In µs
    radicand = 16*R*R-8*R*p+p*p+16*m1Target/(dt*dt*T.Gmax);
    optF     = ceil((-4*R-p+sqrt(radicand))/8); % Solve the m1tmp equation for F
    if strcmpi(func.gradAmp,'scale')
        optF = max(1,optF);
    end
    T.Gmax   = m1Target/(dt*dt*(R + optF)*(p + 4*optF));
    
elseif strcmpi(func.vsType,'BIR8')
    % Due to symmetry can just calculate for one pair of the graients
    m1       = m1Target/2;
    p        = 2*(T.vspad1+T.vspad2+T.RFr+T.r)*1e3/dt; % In µs
    radicand = 8*m1 + dt*dt*T.Gmax*p*p;
    optF     = ceil((-(dt*sqrt(T.Gmax)*(p + 4*R)) + sqrt(radicand))/(4*dt*sqrt(T.Gmax))); % Solve the m1tmp equation for F
    if strcmpi(func.gradAmp,'scale')
        optF = max(1,optF);
    end
    T.Gmax   = m1/(dt*dt*(R + optF)*(p + 2*R + 2*optF));
    
elseif strcmpi(func.vsType,'BIR4')
    p = ceil((T.vspad1+T.vspad2+T.RFr)*1e3/dt);
    radicand = 4*m1Target + dt*dt*T.Gmax*(p+R)*(p+R);
    optF = ceil((-(dt*sqrt(T.Gmax)*(p + 3*R)) + sqrt(radicand))/(2*dt*sqrt(T.Gmax)));
    if strcmpi(func.gradAmp,'scale')
        optF = max(1,optF);
    end
    T.Gmax = m1Target/(dt*dt*(R + optF)*(p + 2*R + optF));
    
elseif strcmpi(func.vsType,'FTVSI')
    % Multiplying Vcut by this factor is a hack to match the "1-crossing" to other methods
    if func.bsinc; kluge = 0.945474;
    else;          kluge = 1.240315; end
    A    = pi/(2*T.Vcut*kluge*T.gamrad*T.Nk*T.Gmax);
    B    = (T.vspad1+T.vspad2+T.RFr/2+T.RFe_2)*1e-3;  % in s
    f    = (sqrt(2*A+(B+T.r*1e-3)^2)-B-3*T.r*1e-3)/2; % in s
    optF = ceil(f*1e6/dt); % Gradient raster points
    if strcmpi(func.gradAmp,'scale')
        optF = max(1,optF);
    end
    delta = T.r*1e-3 + optF*dt*1e-6;
    sep = 2*(2*T.r*1e-3 + optF*dt*1e-6 + B);
    T.Gmax = pi/(2*T.Vcut*kluge*T.gamrad*T.Nk*delta*sep);
    
end

if radicand < 0; error('Complex root: flat top time not possible')
elseif optF < 1; error('f < GUP.');
end

% If we get to here, flatTop is set to optF
T.f = RUP_GRD_ms( optF*dt*1e-3 ); % In ms

% Now calculate the remaining module timings
if     strcmpi(func.vsType,'DRHS');  T = gradtimingsDRHS(T);
elseif strcmpi(func.vsType,'DRHT');  T = gradtimingsDRHS(T);
elseif strcmpi(func.vsType,'BIR8');  T = gradtimingsBIR8(T);
elseif strcmpi(func.vsType,'BIR4');  T = gradtimingsBIR4(T);
elseif strcmpi(func.vsType,'FTVSI'); T = gradtimingsFTVSI(T);
end

    function [T] = gradtimingsDRHS(T)
        T.ta1  = RUP_GRD_ms( T.vspad1                  );
        T.td1  = RUP_GRD_ms( T.ta1 + T.r + T.f + T.r   );
        T.RFr1 = RUP_GRD_ms( T.td1 + T.vspad2          );
        T.ta2  = RUP_GRD_ms( T.RFr1 + T.RFr + T.vspad1 );
        T.td2  = RUP_GRD_ms( T.ta2 + T.r + T.f + T.r   );
        T.ta3  = RUP_GRD_ms( T.td2 + T.vspad2 + 2*T.RFe_2 + T.vspad1 );
        T.td3  = RUP_GRD_ms( T.ta3 + T.r + T.f + T.r   );
        T.RFr2 = RUP_GRD_ms( T.td3 + T.vspad2          );
        T.ta4  = RUP_GRD_ms( T.RFr2 + T.RFr + T.vspad1 );
        T.td4  = RUP_GRD_ms( T.ta4 + T.r + T.f + T.r   );
        T.RFrpad = RUP_GRD_ms( 2*(T.RFe_2+T.vspad1+2*T.r+T.f+T.vspad2) ); % Inter-HS padding
        if abs(T.RFrpad-RUP_GRD_ms(T.RFr2-(T.RFr1+T.RFr))) > T.RFUP; error('RFrpad is wrong'); end
        T.sep23 = RUP_GRD_ms( T.ta3 - T.ta2    ); % Seperation between G2 and G3
        T.RFe2  = RUP_GRD_ms( T.td4 + T.vspad2 ); % Start time of flip up RF
        T.Tvs   = RUP_GRD_ms( T.RFe2 + 2*T.RFe ); % Total VS module time
    end

    function [T] = gradtimingsBIR8(T)
        T.ta1  = RUP_GRD_ms( T.vspad1                  );
        T.td1  = RUP_GRD_ms( T.ta1  + T.r + T.f + T.r  );
        T.RFr1 = RUP_GRD_ms( T.td1  + T.vspad2         );
        T.ta2  = RUP_GRD_ms( T.RFr1 + T.RFr + T.vspad1 );
        T.td2  = RUP_GRD_ms( T.ta2  + T.r + T.f + T.r  );
        T.RFr2 = RUP_GRD_ms( T.td2  + T.vspad2         );
        T.ta3  = RUP_GRD_ms( T.RFr2 + T.RFr + T.vspad1 );
        T.td3  = RUP_GRD_ms( T.ta3  + T.r + T.f + T.r  );
        T.RFr3 = RUP_GRD_ms( T.td3  + T.vspad2         );
        T.ta4  = RUP_GRD_ms( T.RFr3 + T.RFr + T.vspad1 );
        T.td4  = RUP_GRD_ms( T.ta4  + T.r + T.f + T.r  );
        T.RFe2 = RUP_GRD_ms( T.td4 + T.vspad2          ); % Start time of flip up RF
        T.Tvs  = RUP_GRD_ms( T.RFe2 + 2*T.RFe          ); % Total VS module time
    end

    function [T] = gradtimingsBIR4(T)
        T.ta1  = RUP_GRD_ms( T.vspad1                  );
        T.td1  = RUP_GRD_ms( T.ta1  + T.r + T.f + T.r  );
        T.RFr1 = RUP_GRD_ms( T.td1  + T.vspad2         );
        T.ta2  = RUP_GRD_ms( T.RFr1 + T.RFr + T.vspad1 );
        T.td2  = RUP_GRD_ms( T.ta2  + T.r + T.f + T.r  );
        T.RFe2 = RUP_GRD_ms( T.td2 + T.vspad2          ); % Start time of flip up RF
        T.Tvs  = RUP_GRD_ms( T.RFe2 + 2*T.RFe          ); % Total VS module time
    end

    function [T] = gradtimingsFTVSI(T)
        T.ta1  = RUP_GRD_ms( T.vspad1                  );
        T.td1  = RUP_GRD_ms( T.ta1 + T.r + T.f + T.r   );
        T.RFr1 = RUP_GRD_ms( T.td1 + T.vspad2          );
        T.ta2  = RUP_GRD_ms( T.RFr1 + T.RFr + T.vspad1 );
        T.td2  = RUP_GRD_ms( T.ta2 + T.r + T.f + T.r   );
        T.ta3  = RUP_GRD_ms( T.td2 + T.vspad2 + 2*T.RFe_2 + T.vspad1 );
        T.td3  = RUP_GRD_ms( T.ta3 + T.r + T.f + T.r   );
        T.RFr2 = RUP_GRD_ms( T.td3 + T.vspad2          );
        T.ta4  = RUP_GRD_ms( T.RFr2 + T.RFr + T.vspad1 );
        T.td4  = RUP_GRD_ms( T.ta4 + T.r + T.f + T.r   );
        T.RFrpad = RUP_GRD_ms( 2*(T.RFe_2+T.vspad1+2*T.r+T.f+T.vspad2) ); % Inter-HS padding
        if abs(T.RFrpad-RUP_GRD_ms(T.RFr2-(T.RFr1+T.RFr))) > T.RFUP; error('RFrpad is wrong'); end
        T.sep23 = RUP_GRD_ms( T.ta3 - T.ta2      ); % Seperation between G2 and G3
        T.RFe2  = RUP_GRD_ms( T.td4 + T.vspad2   ); % Start time of flip up RF
        T.Tvs   = RUP_GRD_ms( 8*T.RFe2 + 9*T.RFe ); % Total VS module time
    end

end
