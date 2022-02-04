function [PSR,kba,T2a,T2b,R1obs,chi2,chi2p,res,resn] = Analysis_Yarnykh_Full_Fit_DoubleAngle(p,varargin)
%% Analysis_Yarnykh_Full_Fit_DoubleAngle - Finds T1 and Runs the Fitting function for the SP Analysis
% Finds T1 using the Double Angle method from Seth's T1/T2 paper, then prepares the data to run
% the fitting function for finding the full fit qMT
%
% Syntax:  [PSR,kba,T2a,T2b,R1obs,chi2,chi2p,res,resn] = Analysis_Yarnykh_Full_Fit_DoubleAngle(p,varargin)
%
% Inputs:
%   p - Structure which needs the following inputs:
%       p.B1 - B1 correction data (fractional)
%       p.B0 - B0 correction data
%       p.T1flip - FA's for the T1 data in degrees (calculates rest of FA in code)
%       p.Ernst - T1 Data (in order from largest FA to smallest FA)
%       p.T1TR - TR from T1 data (in ms)
%       p.pwMT - pulse length of MT pulse (in s)
%       p.MT_flip - Flip angle of MT saturation pulse (in degrees)
%       p.qMTflip - Flip angle of Excitation pulse (in degrees)
%       p.M - MT data
%       p.TR - TR of MT data (in ms)
% 
% Outputs:
%   PSR - PSR of the MT data
%   kba - kba of the MT data (1/s)
%   T2a - T2 of the free pool for the MT data (s)
%   T2b - T2 of the macromolecular pool for the MT data (s)
%   R1obs - R1obs (in 1/s)
%   chi2 - chi2 goodness of fit
%   chi2p - p-value of chi2 goodness of fit
%   res - residuals
%   resn - normalized residual from lsqnonlin function
%


if nargin > 1
    lineshape = varargin{1};
else
    lineshape = 'super-lorentzian';
end

%% Double Angle - R1obs
corrB1 = double(p.B1/100);

alpha1= double(p.T1flip(1));
alpha2= double(p.T1flip(2));
MFA_Data = double(p.Ernst);
T1TR = double(p.T1TR)*1e-3;

T1obs = TwoAngle_T1_Calc(MFA_Data(1),MFA_Data(2),corrB1,alpha1,alpha2,T1TR);

R1obs = 1/T1obs; % s

if R1obs < 0 || isnan(R1obs)
    PSR = 0;
    kba = 0;
    T2a = 0;
    T2b = 0;
    chi2 = 0;
    chi2p = 0;
    res = 0;
    resn = 0;
    return
end

%% qMT Data Prep

pwMT = double(p.pwMT);
MT_flip = double(p.MT_flip);
qMTflip = double(p.qMTflip);
deltaMT = double(p.deltaMT);
M = double(p.M);
TR = double(p.TR);


gamma = 42.58*2*pi; % Larmor - rad/s-uT


[B1MT,tMT] = philipsRFpulse_FA(MT_flip,pwMT,'am_sg_100_100_0');
B1eMT = CWEqMTPulse(B1MT*corrB1,tMT,pwMT);
thetaEX = ([qMTflip qMTflip]*pi/180)*corrB1;

ts = 1e-3; % s
corrB0 = double(p.B0);
TR = [TR*10^-3, TR*10^-3];

%% 4 Parm Fit

%    PSR   kba    T2a      T2b - Initial Guesses
p0 = [0.15  10  10e-3      10e-6];

[x,chi2,chi2p,res,resn] = fit_SSPulseMT_yarnykh_Full_Fit(p0,M,pwMT,ts,TR,R1obs,thetaEX,B1eMT,deltaMT+corrB0,lineshape);

PSR = x(1);
kba = x(2);
T2a = x(3);
T2b = x(4);

