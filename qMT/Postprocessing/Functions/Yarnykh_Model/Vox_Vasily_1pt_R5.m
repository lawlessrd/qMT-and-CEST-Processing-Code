function Vox_Vasily_1pt_R5(Workspace, deltaMT,pwMT,MT_flip,kba,T2aR1a,T2b,varargin)
%% Vox_Vasily_1pt_R5 - Performs a voxelwise analysis of the Yarnykh method
% This runs the Analysis_Yarnykh_1pt method from the output from SC_Coreg 
% over a large volume in a voxelwise manner.  It will then save out a mat 
% file of the resulting data.
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

[xmax,ymax,zmax,~] = size(qMT_Reg_crop);

if nargin > 7
    PSR_mask = load_nii(varargin{1});
    PSR_mask = PSR_mask.img > 1;
else
    PSR_mask = ones(xmax,ymax,zmax);
end

Parms.T1flip = T1i.info.imgdef.image_flip_angle_in_degrees.uniq;
Parms.MFA = str2double(T1i.info.pardef.Max_number_of_dynamics);
Parms.T1TR = str2double(T1i.info.pardef.Repetition_time_ms);
Parms.qMTflip = qMTi.info.imgdef.image_flip_angle_in_degrees.uniq;
Parms.deltaMT = deltaMT;
Parms.pwMT = pwMT;
Parms.MT_flip = MT_flip;
Parms.kba = kba;
Parms.T2aR1a = T2aR1a;
Parms.T2b = T2b;

Parms.TR = str2double(qMTi.info.pardef.Repetition_time_ms);

PSR_1pt = zeros(xmax,ymax,zmax-2,1);
R1obs = PSR_1pt;
tic;
% Running Middle 10 Slices to be rid of edge effects/partial volume effects
% in Volume
try
    
    for zz = 2:zmax-1
        
        display(['Running Slice ' num2str(zz) ' of ' num2str(zmax)]);
        for ii = 1:xmax
            for jj = 1:ymax
                if PSR_mask(ii,jj,zz) ~= 0
                    
                    display(['At Position ' num2str(ii) ', ' num2str(jj)...
                        ' of Slice ' num2str(zz)]);
                    
                    Parms.MT_flip = MT_flip(1);
                    Parms.deltaMT = deltaMT(1);
                    
                    Parms.B1 = B1_Reg_crop(ii,jj,zz,1);
                    for xx = 1:size(T1_Reg_crop,4);
                        Parms.Ernst(xx) = mean2(T1_Reg_crop(ii,jj,zz,xx));
                    end
                    Parms.B0 = B0_Reg_crop(ii,jj,zz,1);
                    
                    Mz = squeeze(qMT_Reg_crop(ii,jj,zz,:));
                    
                    M = Mz./Mz(end,1);
                    
                    Parms.M = M(1);
                    
                    [PSR_1pt(ii,jj,zz-1),R1obs(ii,jj,zz-1)] = Analysis_Yarnykh_1pt(Parms);
                end
            end
        end
    end
    toc;
    save([Workspace '_with_1pt_Maps']);
    
catch
    toc
    save(sprintf('%s_1pt_Maps_to_slice_%i_of_%i',Workspace,zz,zmax));
end