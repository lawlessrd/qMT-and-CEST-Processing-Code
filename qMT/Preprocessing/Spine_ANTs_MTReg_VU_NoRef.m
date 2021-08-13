function Spine_ANTs_MTReg_VU_NoRef(ScanNumber,Ref,MT,B1,B0,MFA,MT_UpSample,varargin)
%% SPINE_ANTS_MTREG - Runs Full ANTs Registration for Spine MT Data
% Runs Reg_Prep and Reg_ANTs much easier than by using a different script
% for every volume.  Will make producing registered SC data much faster
% than with a script.
%
% Syntax:  Spine_ANTs_MTReg(ScanNumber,Ref,MT,B1,B0,MFA,MT_UpSample,[ImagePath])
%
% Inputs:
%    ScanNumber - The Scan Number Associated with the Data, used for
%    identifying where the data came from, and for saving the Preregistered
%    and Post-Registered Data into a folder.
% 
%    Ref - Filename (in PAR/REC format) of reference Image for Data to be 
%    registered to (mFFE)
% 
%    MT - Filename (in PAR/REC format) of MT Data.
% 
%    B1 - Filename (in PAR/REC format) of B1 Data.
% 
%    B0 - Filename (in PAR/REC format) of B0 Data.
% 
%    MFA - Filename (in PAR/REC format) of MFA Data.
% 
%    MT_UpSample - Flag that specifies if the MT data should be upsampled
%    to mFFE.
%           1 - MT data upsampled to mFFE (and all other data resized to
%               mFFE)
%           0 - Ref sampled to MT data (and all other data resized to MT)
% 
%    ImagePath - (optional) Path to images if PAR/REC files are not located in same
%    directory as current folder.
%            Default - ImagePath = pwd
%    
%
% Example: 
% ImagePath = '/Volumes/ms_sas/3T/Spine/controls/Smith_215347_20150508';
% ScanNumber = '215347_FF';
% 
% MT = 'Smith_215347_WIP_PulseMT_16_dyn_SENSE_10_1.PAR';
% Ref = 'Smith_215347_sWIP_mFFE_0.65_SENSE_4_3.PAR';
% B1 = 'Smith_215347_WIP_WIP_B1_Map_Yarnykh_CLEAR_5_1.PAR';
% B0 = 'Smith_215347_WIP_B0_Map_In-phase_6_1.PAR';
% MFA = 'Smith_215347_WIP_MFA_T1_SENSE_9_1.PAR';
%
% 
% Spine_ANTs_MTReg(ScanNumber,Ref,MT,B1,B0,MFA,1,ImagePath)
%
% Subfunctions: Reg_Prep, Reg_ANTs, FSL,ANTs
%
% Author:  smithak4
% Date:    17-Feb-2016
% Version: 1.0
% Changelog:
%
% 17-Feb-2016 - initial creation
%
%------------- BEGIN CODE --------------


%% Define Inputs
if nargin > 7
    ImagePath = varargin{1};
else
    ImagePath = pwd;
end

%% Prep Volumes for Registration by Saving them out to NIFTI

MFA_detrend = [1.8287    1.6001    1.2828    1.1001    1.0000 1.1725];

PrepTime = tic;
% 
if MT_UpSample
    
    Reg_Prep_vu_B1nii(Ref,Ref,ScanNumber,'mFFE',ImagePath);
    Reg_Prep_vu_B1nii(B1,Ref,ScanNumber,'B1',ImagePath);
    Reg_Prep_vu_B1nii(MT,Ref,ScanNumber,'MT',ImagePath);
    Reg_Prep_vu_B1nii(MFA,Ref,ScanNumber,'MFA',ImagePath,[],MFA_detrend);
    Reg_Prep_vu_B1nii(B0,Ref,ScanNumber,'B0',ImagePath);
    
else
    
%     MT_info = loadPARREC(sprintf('%s/%s',ImagePath,MT));
%     f= MT_info.imgdef.recon_resolution_x_y.uniq;
    resx = 256;%f(1);
    resy = 256;%f(2);

    Reg_Prep_vu_B1nii(Ref,Ref,ScanNumber,'mFFE',ImagePath,[resx, resy]);
    Reg_Prep_vu_B1nii(B1,Ref,ScanNumber,'B1',ImagePath,[resx, resy]);
    Reg_Prep_vu_B1nii(MT,Ref,ScanNumber,'MT',ImagePath,[resx, resy]);
    Reg_Prep_vu_B1nii(MFA,Ref,ScanNumber,'MFA',ImagePath,[resx, resy],MFA_detrend);
    Reg_Prep_vu_B1nii(B0,Ref,ScanNumber,'B0',ImagePath,[resx, resy]);

end


fprintf('\nTime to Perform Registration Prep is: %3.4f minutes \n', toc(PrepTime)/60);


%% Register all Volumes to Reference

AntsTime = tic;

Reg_ANTs_vu_NoRef(ScanNumber,'MT');
Reg_ANTs_vu_NoRef(ScanNumber,'B0');
Reg_ANTs_vu_NoRef(ScanNumber,'B1');
Reg_ANTs_vu_NoRef(ScanNumber,'MFA');

fprintf('\nTime to Perform Registration is: %3.4f minutes \n', toc(AntsTime)/60);

%% Load Volumes and Perform Final Crop

% % Check if scanID is a character
% if ~ischar(ScanNumber)
%     ScanNumber = num2str(ScanNumber);
% end
% 
% unix(sprintf('gunzip %s_Registration/%s*.gz',ScanNumber,ScanNumber));
% 
% tmp = load_nii(sprintf('%s_Registration/%s_MT_Reg.nii',ScanNumber,ScanNumber));
% qMT_Reg = tmp.img;
% % qMTi.info = loadPARREC(sprintf('%s/%s',ImagePath,MT));
% 
% tmp = load_nii(sprintf('%s_Registration/%s_MFA_Reg_to_MFA.nii',ScanNumber,ScanNumber));
% T1_Reg = tmp.img;
% % T1i.info = loadPARREC(sprintf('%s/%s',ImagePath,MFA));
% 
% tmp = load_nii(sprintf('%s_Registration/%s_B1_Reg.nii',ScanNumber,ScanNumber));
% B1_Reg = tmp.img;
% % B1i.info = loadPARREC(sprintf('%s/%s',ImagePath,B1));
% 
% tmp = load_nii(sprintf('%s_Registration/%s_B0_Reg.nii',ScanNumber,ScanNumber));
% B0_Reg = tmp.img;
% % B0i.info = loadPARREC(sprintf('%s/%s',ImagePath,B0));
% 
% tmp = load_nii(sprintf('%s_Registration/%s_Reference.nii',ScanNumber,ScanNumber));
% Ref_Reg = tmp.img;
% % Refi.info = loadPARREC(sprintf('%s/%s',ImagePath,Ref));
% 
% unix(sprintf('gzip %s_Registration/%s*.nii',ScanNumber,ScanNumber));
% 
% [rmax,~,zmax] = size(Ref_Reg);
% 
% rhalf = round(rmax/2);
% 
% Crop = floor(rmax*.35);
% 
% Ref_Reg_crop = zeros(2*Crop,2*Crop,zmax);
% B1_Reg_crop = zeros(2*Crop,2*Crop,zmax);
% B0_Reg_crop = zeros(2*Crop,2*Crop,zmax);
% 
% qMT_ndyn = size(qMT_Reg,4);
% MFA_ndyn = size(T1_Reg,4);
% 
% qMT_Reg_crop = zeros(2*Crop,2*Crop,zmax,qMT_ndyn);
% T1_Reg_crop = zeros(2*Crop,2*Crop,zmax,MFA_ndyn);
% 
% for kk = 1:zmax
%     Ref_Reg_crop(:,:,kk)=imcrop(Ref_Reg(:,:,kk),[rhalf-Crop rhalf-Crop 2*Crop-1 2*Crop-1]);
%     B1_Reg_crop(:,:,kk)=imcrop(B1_Reg(:,:,kk),[rhalf-Crop rhalf-Crop 2*Crop-1 2*Crop-1]);
%     B0_Reg_crop(:,:,kk)=imcrop(B0_Reg(:,:,kk),[rhalf-Crop rhalf-Crop 2*Crop-1 2*Crop-1]);
%     
%     for tt1 = 1:qMT_ndyn
%         qMT_Reg_crop(:,:,kk,tt1)=imcrop(qMT_Reg(:,:,kk,tt1),[rhalf-Crop rhalf-Crop 2*Crop-1 2*Crop-1]); 
%     end
%     for tt2 = 1:MFA_ndyn
%         T1_Reg_crop(:,:,kk,tt2)=imcrop(T1_Reg(:,:,kk,tt2),[rhalf-Crop rhalf-Crop 2*Crop-1 2*Crop-1]);
%     end
% end
% 
% clear kk tt1 tt2 tmp rhalf ans
% 
% 
% save(sprintf('MT_CoReg_%s',ScanNumber));

fprintf('\n--------------------------------------- \n');
fprintf('   Finished Registration of %s. \n',ScanNumber);
fprintf('--------------------------------------- \n');


%------------- END OF CODE --------------