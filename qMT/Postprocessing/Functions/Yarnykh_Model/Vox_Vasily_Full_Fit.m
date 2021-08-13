function Vox_Vasily_Full_Fit(Workspace, deltaMT,pwMT,MT_flip,varargin)
%% Vox_Vasily_Full_Fit - Performs a voxelwise analysis of the Yarnykh method
% This runs the Analysis_Yarnykh_Full_Fit method from the output from SC_Coreg 
% over a large volume in a voxelwise manner.  It will then save out a mat 
% file of the resulting data.
%
% Syntax:  Vox_Vasily_Full_Fit(Workspace, deltaMT,pwMT,MT_flip,[Mask])
%
% Inputs:
%   Workspace - The mat file created by SC_Coreg which will then be used to
%       perform the analysis.  This will also be used to name the file that
%       is output at the end of the analysis.
%   deltaMT - The offsets (in Hz) for the MT data
%   pwMT - Pulse width (in s) of the MT pulse
%   MT_flip - Pulse power (in degrees) of the MT pulse
% 
%   Optional Inputs:
%       Mask - Mask of the volume of interest.
% 
%

load(Workspace);


[xmax,ymax,zmax,tmax] = size(qMT_Reg_crop);


if nargin > 4
    mask = varargin{1};
else
    mask = ones(xmax,ymax,zmax);
end

Parms.T1flip = mean(T1i.Parms.tags(:,36));
Parms.MFA = T1i.Parms.max_dynamics;
Parms.T1TR = T1i.Parms.repetition_time; %ms
Parms.qMTflip = mean(qMTi.Parms.tags(:,36));
Parms.deltaMT = deltaMT;
Parms.pwMT = pwMT;
Parms.MT_flip = MT_flip;

scan = length(Parms.deltaMT);
Parms.TR = qMTi.Parms.repetition_time;

MTR = zeros(xmax,ymax,zmax-2,tmax);
PSR = zeros(xmax,ymax,zmax-2);
kba = PSR;
T2a = PSR;
T2b = PSR;
R1obs = PSR;
chi2 = PSR;
chi2p = PSR;
tic;

% Running Middle Slices to be rid of edge effects/partial volume effects
% in Volume

for kk = 2:zmax-1
    display(['Running Slice ' num2str(kk) ' of ' num2str(zmax-2)]);
    
    for ii = 1:xmax
        for jj = 1:ymax
            if mask(ii,jj,kk) ~= 0
                display(['At Position ' num2str(ii) ', ' num2str(jj)...
                    ' of Slice ' num2str(kk)]);
                
                Parms.B1 = B1_Reg_crop(ii,jj,kk,1);
                for xx = 1:size(T1_Reg_crop,4);
                    Parms.Ernst(xx) = mean2(T1_Reg_crop(ii,jj,kk,xx));
                end
                Parms.B0 = B0_Reg_crop(ii,jj,kk,1);
                
                
                for xx = 1:size(qMT_Reg_crop,4)
                    Mz(xx,1) = qMT_Reg_crop(ii,jj,kk,xx);
                end
                
                M(1:scan,1) = Mz(1:scan);
                M(1:scan,2) = Mz(scan+1:end);
                
                Parms.M(:,1) = M(:,1)./M(end,1);
                Parms.M(:,2) = M(:,2)./M(end,2);
                
                [PSR(ii,jj,kk-1),...
                    kba(ii,jj,kk-1),...
                    T2a(ii,jj,kk-1),...
                    T2b(ii,jj,kk-1),...
                    R1obs(ii,jj,kk-1),...
                    chi2(ii,jj,kk-1),...
                    chi2p(ii,jj,kk-1),] = Analysis_Yarnykh_Full_Fit(Parms);
                
                for tt = 1:tmax
                    % Solving for each offset independently
                    if tt <= 8
                        MTR(ii,jj,kk-1,tt) = 1-qMT_Reg_crop(ii,jj,kk,tt)/qMT_Reg_crop(ii,jj,kk,8);
                    else
                        MTR(ii,jj,kk-1,tt) = 1-qMT_Reg_crop(ii,jj,kk,tt)/qMT_Reg_crop(ii,jj,kk,tmax);
                    end
                    
                    % To correct for Offsets less than the Ref Being over 1
                    if MTR(ii,jj,kk-1,tt) < 0
                        MTR(ii,jj,kk-1,tt) = 0;
                    end
                    
                end
            end
            
        end
    end
end
time = toc;

disp(['The time to find a full fit map for the data is: ' num2str(time)]);
save([Workspace '_with_PSR_Maps']);
