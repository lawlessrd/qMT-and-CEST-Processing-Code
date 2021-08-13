function Vox_Vasily_1pt(Workspace, deltaMT,pwMT,MT_flip,kba,T2aR1a,T2b,varargin)
%% Vox_Vasily_1pt - Performs a voxelwise analysis of the Yarnykh method
% This runs the Analysis_Yarnykh_Full_Fit method from the output from SC_Coreg 
% over a large volume in a voxelwise manner.  It will then save out a mat 
% file of the resulting data.
% 
% This has syntax built for data collected before R5.  (legacy to be able
% to run old data).
%
% Syntax:  Vox_Vasily_1pt_R5(Workspace, deltaMT,pwMT,MT_flip,kba,T2aR1a,T2b,[Mask])
%
% Inputs:
%   Workspace - The mat file created by SC_Coreg which will then be used to
%       perform the analysis.  This will also be used to name the file that
%       is output at the end of the analysis.
%   deltaMT - The offsets (in Hz) for the MT data
%   pwMT - Pulse width (in s) of the MT pulse
%   MT_flip - Pulse power (in degrees) of the MT pulse
%   kba - exchange rate from full fit analysis (in 1/s).
%   T2aR1a - T2aR1a from full fit analysis.
%   T2b - T2b from full fit analysis (in s).
% 
%   Optional Inputs:
%       Mask - Mask of the volume of interest.
% 
%

load(Workspace);

[xmax,ymax,zmax,tmax] = size(qMT_Reg_crop);

if nargin > 7
    PSR_mask = varargin{1};
else
    PSR_mask = ones(xmax,ymax,zmax);
end

Parms.T1flip = mean(T1i.Parms.tags(:,36));
Parms.MFA = T1i.Parms.max_dynamics;
Parms.T1TR = T1i.Parms.repetition_time; %ms
Parms.qMTflip = mean(qMTi.Parms.tags(:,36));
Parms.deltaMT = deltaMT;
Parms.pwMT = pwMT;
Parms.MT_flip = MT_flip;

Parms.kba = kba;
Parms.T2aR1a = T2aR1a;
Parms.T2b = T2b;

scan = length(Parms.deltaMT);
Parms.TR = qMTi.Parms.repetition_time;

numoffsets = length(MT_flip)*length(deltaMT);

PSR_1pt = zeros(xmax,ymax,zmax-2,numoffsets);
R1obs_1pt = zeros(xmax,ymax,zmax-2);

tic;
for kk = 1:numoffsets
    display(['Running Offset ' num2str(kk) ' of ' num2str(numoffsets)]);
    
    % Running Middle 10 Slices to be rid of edge effects/partial volume effects
    % in Volume
    
    for zz = 2:zmax-1
        
        display(['Running Slice ' num2str(zz) ' of ' num2str(zmax)]);
        for ii = 1:xmax
            for jj = 1:ymax
                if PSR_mask(ii,jj,zz) ~= 0
                    
                    display(['At Position ' num2str(ii) ', ' num2str(jj)...
                        ' of Slice ' num2str(zz)]);
                    
                    off_length = length(deltaMT);
                    if kk <= off_length
                        Parms.MT_flip = MT_flip(1);
                        Parms.deltaMT = deltaMT(kk);
                    else
                        Parms.MT_flip = MT_flip(2);
                        Parms.deltaMT = deltaMT(kk-off_length);
                    end
                    
                    Parms.B1 = B1_Reg_crop(ii,jj,zz,1);
                    for xx = 1:size(T1_Reg_crop,4);
                        Parms.Ernst(xx) = mean2(T1_Reg_crop(ii,jj,zz,xx));
                    end
                    Parms.B0 = B0_Reg_crop(ii,jj,zz,1);
                    
                    
                    for xx = 1:size(qMT_Reg_crop,4)
                        Mz(xx,1) = qMT_Reg_crop(ii,jj,zz,xx);
                    end
                    
                    M(1:scan,1) = Mz(1:scan);
                    
                    if length(MT_flip) ==2
                        M(1:scan,2) = Mz(scan+1:end);
                    end
                    
                    M2(:,1) = M(:,1)./M(end,1);
                    
                    if length(MT_flip) ==2
                        M2(:,2) = M(:,2)./M(end,2);
                    end
                    
                    
                    
                    if kk <= off_length
                        Parms.M = M2(kk,1);
                    else
                        Parms.M = M2(kk-off_length,2);
                    end
                    
                    [PSR_1pt(ii,jj,zz-1,kk),R1obs_1pt(ii,jj,zz-1,kk)] = Analysis_Yarnykh_1pt(Parms);
                else
                    PSR_1pt(ii,jj,zz-1,kk) = 0;
                    R1obs_1pt(ii,jj,zz-1,kk) = 0;
                end
            end
        end
    end
end
time = toc;

save([Workspace '_with_1pt_Maps']);