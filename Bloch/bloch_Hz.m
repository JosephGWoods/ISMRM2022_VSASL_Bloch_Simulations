%   Bloch simulator that works in units of Hz (to treat all nuclei equally).
%
%   [mx,my,mz] = bloch_Hz(b1,gr,tp,t1,t2,df,dp,dv,mode,mx,my,mz,spoiled,offset)
%
%   Bloch simulation of rotations due to B1, gradient and off-resonance,
%   including relaxation effects. At each time point, the rotation matrix
%   and decay matrix are calculated. Simulation can simulate the
%   steady-state if the sequence is applied repeatedly, or the
%   magnetization starting at m0.
%
%   in:
%       b1   - (Mx1) RF pulse (Hz).  Can be complex.
%       gr   - (Mx1,2,3) 1,2 or 3-dimensions are linear gradients (Hz/cm) x,y,z.
%       tp   - (Mx1) time duration of each b1 and gr point (seconds),
%              or 1x1 time step if constant for all points
%              or monotonically INCREASING endtime of each interval.
%       t1   - T1 relaxation time (seconds).
%       t2   - T2 relaxation time (seconds).
%       df   - (Fx1) Array of off-resonance frequencies (Hz).
%       dp   - (Px1,2,or 3) Array of spatial positions (cm).
%       dv   - (Vx1,2,or 3) Array of velocities (cm/s)
%
%	(optional)
%       mode	 - Bitmask mode:
%                  Bit 0:  0-Simulate from start or M0, 1-Steady State
%                  Bit 1:  1-Record m at time points.  0-just end time.
%                    mode = 0 -> Bit0=0, Bit1=0
%                    mode = 1 -> Bit0=1, Bit1=0
%                    mode = 2 -> Bit0=0, Bit1=1
%                    mode = 3 -> Bit0=1, Bit1=1
%       mx,my,mz - (PxFxV) arrays of starting magnetization, where N
%                  is the number of frequencies, P is the number
%                  of spatial positions, and V is the number of velocities.
%
%   out:
%       mx,my,mz - (PxFxV) arrays of the resulting magnetization components
%                  at each position and frequency.
%
%   B. Hargreaves   Nov 2003. (Downloaded in Feb 2013.)
%   M. Robson       REVERSED SIGN OF GYROMAGNETIC RATIO.
%   C. Rodgers      Feb 2013. Add flag to control verbosity.
%                   Don't crash Matlab if called with insufficient
%                   parameters.
%                   SWAPPED SIGN OF GAMMA to match updated code on web.
%                   Jun 2013. Convert from G --> Hz units.
%   J.G. Woods      Apr 2020. Added support for simulating flow with constant
%                   velocities.
%                   Nov 2020. Don't crash Matlab if required parameters are
%                   empty.
%
% SIGN convention. See M. Levitt. "Basics of Nuclear Magnetic Resonance". Page 250.
% An 90-x pulse (i.e. RF with 0 phase) moves spins from the z-axis to the
% MINUS y-axis. i.e. z --> -y.
%
% The updated code here adheres to that convention.

function varargout = bloch_Hz(varargin)
warning('This code is a MEX file that must be compiled before use. Attempting to do that now...');

oldDir = pwd();
c = onCleanup(@() cd(oldDir));

try
    myDir  = fileparts(mfilename('fullpath'));
    myFunc = mfilename();
    myFile = sprintf('%s.c',myFunc);
    
    cd(myDir)
    mex -setup C;
    mex(myFile,'-compatibleArrayDims')
    
catch ME
    error('Error automatically compiling MEX. Must be done by hand.')
end

clear c; % Return to previous folder
rehash;  % Refresh function and file system path caches
warning('Compilation done. Continuing...')

% Run compiled bloch_Hz
myHandle = str2func(myFunc);
[varargout{1:nargout}] = myHandle(varargin{:});
