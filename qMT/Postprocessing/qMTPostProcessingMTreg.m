%% qMT Post Processing 
% RDL 
%% Set paths

% FSL
setenv('FSLOUTPUTTYPE','NIFTI_GZ')
Path = getenv('PATH');
if ~contains(Path,'Users/lawlesrd/Documents/fsl/bin')
    setenv('PATH', [Path ':/Users/lawlesrd/Documents/fsl/bin']);  %Edit rdl 4/19/2018 getenv was just PATH
     % setenv('PATH', [Path ':/usr/local/fsl/bin']);  %home version
end

% ANTS
if ~contains(Path,'/Users/lawlesrd/install/bin/')
    setenv('PATH', [Path ':/Users/lawlesrd/install/bin/']);  %Edit rdl 4/19/2018 getenv was just PATH
      %setenv('PATH', [Path ':/opt/ants/bin/']);  %home version
end
clear PATH

% SCT
if ~contains(Path,'/Users/dylanlawless/sct_5.3.0/bin')
    setenv('PATH', [Path ':/Users/dylanlawless/sct_5.3.0/bin']);  %Edit rdl 4/19/2018 getenv was just PATH
      %setenv('PATH', [Path ':/opt/ants/bin/']);  %home version
end
clear PATH


antsPath = '/Users/lawlesrd/install/bin/';

%% Import scans

tic;
home=pwd;
out = regexp(home, '\d+', 'match');
ScanNumber = [out{1}];
imgPath = sprintf('%s/%s_Registration',home,ScanNumber);
cd(imgPath);

% Note: Make sure we are using magnitude B1

a=dir(sprintf('%s_MT_Reg.nii.gz',ScanNumber));
b=dir(sprintf('%s_MFA_Reg_to_MFA.nii.gz',ScanNumber));
c=dir(sprintf('%s_B0_Reg.nii.gz',ScanNumber));
d=dir(sprintf('%s_B1_Reg_to_MT.nii.gz',ScanNumber));
e=dir(sprintf('%s_Reference.nii.gz',ScanNumber)); % mFFE

g=dir(sprintf('%s_B1_Reg_to_MFA.nii.gz',ScanNumber));

MT_file=a.name;
MFAtoMFA_file=b.name;
B0_file=c.name;
B1toMT_file=d.name;
Ref_file=e.name;

B1toMFA_file=g.name;

qMT=niftiread(MT_file);
B0=niftiread(B0_file);
B1toMT=niftiread(B1toMT_file);
B1toMFA=niftiread(B1toMFA_file);
T1toT1=niftiread(MFAtoMFA_file);

Ref=niftiread(Ref_file);


%% Input parameters
% %MT GORE
% p.TR = 42.3; %MT TR (ms)
% p.qMTflip = 6; %Readout MT FA (degrees)
% p.MT_flip = [360 820]; %MT prepulse (MT power) (degrees) 
% p.pwMT = [20e-3 20e-3]; %Prepulse duration (s)
% p.deltaMT = [1000 1500 2000 2500 8000 16000 32000 100000]; %MT offsets (Hz)
% 
% %T1 MFA
% p.T1flip = 30; % MFA flip angle (degrees)
% p.MFA = 6; % # of MFA dynamics (FAs)
% p.T1TR = 20; % T1 TR (ms)

%MT
p.TR = 50; %MT TR (ms)
p.qMTflip = 6; %Readout MT FA (degrees)
p.MT_flip = [360 820]; %MT prepulse (MT power) (degrees) 
p.pwMT = [20e-3 20e-3]; %Prepulse duration (s)
p.deltaMT = [1000 1500 2000 2500 8000 16000 32000 100000]; %MT offsets (Hz)

%T1 MFA
p.T1flip = 30; % MFA flip angle (degrees)
p.MFA = 6; % # of MFA dynamics (FAs)
p.T1TR = 50; % T1 TR (ms)

%UTSW
% p.TR = 57; %MT TR (ms)
% p.qMTflip = 6; %Readout MT FA (degrees)
% p.MT_flip = [600 1000]; %MT prepulse (MT power) (degrees) 
% p.pwMT = [20e-3 20e-3]; %Prepulse duration (s)
% p.deltaMT = [1000 1500 2000 2500 8000 16000 32000 100000]; %MT offsets (Hz)
% 
% %T1 MFA
% p.T1flip = 30; % MFA flipfigure() angle (degrees)
% p.MFA = 6; % # of MFA dynamics (FAs)
% p.T1TR = 20; % T1 TR (ms)

%% Normalization
% Split qMT images by power and normalize to 100 kHz offset

M1 = qMT(:,:,:,1:8);
M2 = qMT(:,:,:,9:end);

for jj = 1 : size(M1,3)
    for ii = 1 : size(M1,4)
        M1(:,:,jj,ii) = M1(:,:,jj,ii)./M1(:,:,jj,end);
        M2(:,:,jj,ii) = M2(:,:,jj,ii)./M2(:,:,jj,end);
    end
end

M1(isinf(M1)) = 0 ;
M1(isnan(M1)) = 0 ;

M2(isinf(M2)) = 0 ;
M2(isnan(M2)) = 0 ;

%% FITTING
R1obs = zeros(size(B0));
PSR = zeros(size(B0));
R1obsPrereg = zeros(size(B0));
kba = zeros(size(B0));
T2a = zeros(size(B0));
T2b = zeros(size(B0));
T1a = zeros(size(B0));
chi2 = zeros(size(B0));
chi2p = zeros(size(B0));
resn = zeros(size(B0));
res = 0;

Mn = zeros(2,size(M1,4));

%% Create mask or segment
% Creates a simple cylindrical mask using Spinal Cord Toolbox

% Note: might be best to process entire image

% Set mask size
maskSize = 50;

unix(sprintf('sct_create_mask -i %s_MT_Reg.nii.gz -p center -size %0.0f -f cylinder -o %s_mask_%0.0f.nii.gz'...
    ,ScanNumber,maskSize,ScanNumber,maskSize));

mask = niftiread(sprintf('%s_mask_%0.0f.nii.gz',ScanNumber,maskSize));

%% Calculate T1 maps
% Requirements: MFA registered to first FA, B1 registered to first FA of MFA

for s = 1:size(M1,3)
    
    [row, col] = find(mask(:,:,s));
    
    fprintf('Slice %g: \n', s)
    fprintf('Number of voxels: %g\n',length(row));
    
    for ii = drange(1:length(row))
        
        if mod(ii/1000,1) == 0
            fprintf('%g/%g: \n',ii, length(row));
        end
        % MFA - R1obs - (Ke's Code)
        p.B1toMT = squeeze(B1toMT(row(ii), col(ii), s));
        p.B1toT1 = squeeze(B1toMFA(row(ii), col(ii), s));
        p.Ernst2MFA = squeeze(T1toT1(row(ii), col(ii),s,:))';
        
        
        corrB1toMFA = double(p.B1toT1/100);
        
        T1flip = double(p.T1flip);
        MFA = double(p.MFA);
        MFAtoMFA_Data = double(p.Ernst2MFA);
        T1TR = double(p.T1TR);
        
        T1toT1obs = T1_Calc(T1flip,MFA,MFAtoMFA_Data,corrB1toMFA,T1TR);

        R1obsPrereg(row(ii), col(ii),s)= 1/T1toT1obs*1000;

    end
end


%% Register R1map to MT using FLIRT
% 
 MFAnii = niftiinfo(sprintf('%s_Prereg/%s_MFA_vol_1.nii.gz',ScanNumber,ScanNumber));
% 
 niftiwrite(R1obsPrereg,'R1obs.nii');
% 
 in_vol = sprintf('R1obs.nii');
 trans_mat = sprintf('%s_Prereg/%s_T1_vol1_to_MT.xfm',ScanNumber,ScanNumber);
 target_vol = sprintf('%s_Prereg/%s_MT_vol_1.nii.gz',ScanNumber,ScanNumber);
out_volR1 = sprintf('%s_R1obs2MT.nii.gz',ScanNumber);

input_trans = sprintf('flirt -ref %s -in %s -2D -applyxfm -init %s -out %s',target_vol,in_vol,trans_mat,out_volR1);

%disp(['Register T1 Data for: ' in_vol ' to ' target_vol])
unix(input_trans);

R1obs = niftiread(out_volR1);

%% MT Processing
for s = 1:size(M1,3)
    
    [row, col] = find(mask(:,:,s));
    
    fprintf('Slice %g: \n', s)
    fprintf('Number of voxels: %g\n',length(row));
    
    for ii = drange(1:length(row))
        
        if mod(ii/1000,1) == 0
            fprintf('%g/%g: \n',ii, length(row));
        end
        
        Mn1 = squeeze(M1(row(ii), col(ii), s, : ))';
        Mn2 = squeeze(M2(row(ii), col(ii), s, : ))';
        Mn(1,:) = Mn1;
        Mn(2,:) = Mn2;
        
        p.B1 = squeeze(B1toMT(row(ii), col(ii), s));
%        p.B1toT1 = squeeze(B1toMFA(row(ii), col(ii), s));
        p.Ernst = squeeze(R1obs(row(ii), col(ii),s,:))';
        p.B0 = squeeze(B0(row(ii), col(ii), s));
        p.M = Mn;
        [PSR(row(ii), col(ii), s),kba(row(ii), col(ii), s),T2a(row(ii),col(ii), s)...
            ,T2b(row(ii), col(ii), s),R1obs(row(ii), col(ii),s),chi2(row(ii), col(ii), s)...
            ,chi2p(row(ii), col(ii), s),res,resn(row(ii), col(ii), s)] = Analysis_Yarnykh_Full_Fit(p);
%         fprintf('PSR=%g, kba=%g, T2a=%g, T2b=%g, T1obs=%g, T1a=%g, chi2=%g, chi2p=%g, resn=%g, sigma=%g \n',PSR(row(ii), col(ii), s),kba(row(ii), col(ii), s),T2a(row(ii), col(ii), s),T2b(row(ii), col(ii), s),T1obs(row(ii), col(ii), s),T1a(row(ii), col(ii), s),chi2(row(ii), col(ii), s),chi2p(row(ii), col(ii), s),resn(row(ii), col(ii), s),sigma(row(ii), col(ii), s));
    end
    
     MPF = PSR ./ (1 + PSR);
     kab = PSR .* kba;
    save('qMT_results_full');
end

%% Calculate MTR and save as nii

    %% Isolate S0 (100 kHz offset), and S_MT (2 and 2.5 kHz offsets)
    % For both powers
    
    S0_low = M1(:,:,:,end);
    SMT_2_low = M1(:,:,:,3);
    SMT_3_low = M1(:,:,:,4);
    
    S0_high = M2(:,:,:,end);
    SMT_2_high = M2(:,:,:,3);
    SMT_3_high = M2(:,:,:,4);
    
    % Calulate MTR using both offsets and both powers
    MTR_2_low = (S0_low - SMT_2_low) ./ S0_low;
    MTR_2_high= (S0_high - SMT_2_high) ./ S0_high;
    
    MTR_3_low = (S0_low - SMT_3_low) ./ S0_low; % best for SC MTR calculations
    MTR_3_high= (S0_high - SMT_3_high) ./ S0_high;
    
    
    
    % Save variables
    
    save(sprintf('%s_MTR',ScanNumber))



%% Move T1obs and PSR to mFFE space
% Steps:
% 1. Resample PSR and R1obs to 14 slices
% 2. Get transform from MT to mFFE
% 3. Apply transform

%% Resample PSR, MT to 14 slices Save results

MTGvol = sprintf('%s_Prereg/%s_MT_Gvol_1.nii.gz',ScanNumber,ScanNumber);
MTvol = sprintf('%s_Prereg/%s_MT_vol_1.nii.gz',ScanNumber,ScanNumber);
Refvol = sprintf('%s_Prereg/%s_mFFE_Gvol_1.nii.gz',ScanNumber,ScanNumber);
T1vol = out_volR1;

MTr = sprintf('%s_Prereg/%s_MT_vol_1_resample.nii.gz',ScanNumber,ScanNumber);
MTGr = sprintf('%s_Prereg/%s_MT_Gvol_1_resample.nii.gz',ScanNumber,ScanNumber);
PSRr = sprintf('%s_Prereg/%s_PSR2MT_resample.nii.gz',ScanNumber,ScanNumber);
R1r = sprintf('%s_Prereg/%s_R12MT_resample.nii.gz',ScanNumber,ScanNumber);
MTRr = sprintf('%s_Prereg/%s_MTR2MT_resample.nii.gz',ScanNumber,ScanNumber);


RefInfo = niftiinfo(Refvol);

% save PSR to nii
PSRnii = niftiinfo(MTvol);
niftiwrite(PSR,sprintf('%s_PSR2MT.nii',ScanNumber),PSRnii);
unix(sprintf('gzip %s_PSR2MT.nii',ScanNumber));

% save MTR to nii
MTRnii = niftiinfo(MTvol);
niftiwrite(MTR_3_low,sprintf('%s_MTR2MT.nii',ScanNumber),MTRnii);
unix(sprintf('gzip %s_MTR2MT.nii',ScanNumber));

% cd(home)
temp_folder = 'ANTs_Reg';
mkdir(temp_folder);

% Resample MT_vol1, MT_Gvol1, PSR2MT, and R1obs2MT to match mFFE
% (150x150x14)
unix(['sct_resample -i ' MTvol ' -vox ' num2str(RefInfo.ImageSize(1)) ...
    'x' num2str(RefInfo.ImageSize(2)) 'x' num2str(RefInfo.ImageSize(3)) ...
    ' -x nn -o ' MTr ]);
unix(['sct_resample -i ' MTGvol ' -vox ' num2str(RefInfo.ImageSize(1)) ...
    'x' num2str(RefInfo.ImageSize(2)) 'x' num2str(RefInfo.ImageSize(3)) ...
    ' -x nn -o ' MTGr ]);
unix(['sct_resample -i ' ScanNumber '_PSR2MT.nii.gz -vox ' num2str(RefInfo.ImageSize(1)) ...
    'x' num2str(RefInfo.ImageSize(2)) 'x' num2str(RefInfo.ImageSize(3)) ...
    ' -x nn -o ' PSRr ]);
unix(['sct_resample -i ' T1vol ' -vox ' num2str(RefInfo.ImageSize(1)) ...
    'x' num2str(RefInfo.ImageSize(2)) 'x' num2str(RefInfo.ImageSize(3)) ...
    ' -x nn -o ' R1r ]);
unix(['sct_resample -i ' ScanNumber '_MTR2MT.nii.gz -vox ' num2str(RefInfo.ImageSize(1)) ...
    'x' num2str(RefInfo.ImageSize(2)) 'x' num2str(RefInfo.ImageSize(3)) ...
    ' -x nn -o ' MTRr ]);

volume =1;

% Break Reference (mFFE) into individual slices

unix(sprintf('fslsplit %s_Prereg/%s_mFFE_vol_1.nii.gz %s/ref_ -z',...
    ScanNumber,ScanNumber,temp_folder));

unix(sprintf('fslsplit %s_Prereg/%s_mFFE_Gvol_1.nii.gz %s/ref_G_ -z',...
    ScanNumber,ScanNumber,temp_folder));


%Break outputs into individual slices
unix(sprintf('fslsplit %s %s/PSRslice_ -z',...
    PSRr,temp_folder));
unix(sprintf('fslsplit %s %s/R1slice_ -z',...
    R1r,temp_folder));
unix(sprintf('fslsplit %s %s/MTRslice_ -z',...
    MTRr,temp_folder));

% Break Image into individual slices
unix(sprintf('fslsplit %s %s/im_vol%i_ -z',...
    MTr,temp_folder,volume));

unix(sprintf('fslsplit %s %s/im_Gvol%i_ -z',...
    MTGr,temp_folder,volume));

num_slices = length(dir(sprintf('%s/ref_0*.nii.gz',temp_folder)));

for ii=1:num_slices
    fprintf('\tProcessing for slice number: %i\n', ii);
    
    if ii-1 < 10
        ref_gauss = sprintf('%s/ref_G_000%i.nii.gz',temp_folder,ii-1);
        ref = sprintf('%s/ref_000%i.nii.gz',temp_folder,ii-1);
        Im_gauss = sprintf('%s/im_Gvol%i_000%i.nii.gz',temp_folder,volume,ii-1);
        Im = sprintf('%s/im_vol%i_000%i.nii.gz',temp_folder,volume,ii-1);
        ImPSR = sprintf('%s/PSRslice_000%i.nii.gz',temp_folder,ii-1);
        ImR1 = sprintf('%s/R1slice_000%i.nii.gz',temp_folder,ii-1);
        ImMTR = sprintf('%s/MTRslice_000%i.nii.gz',temp_folder,ii-1);
        
    else
        ref_gauss = sprintf('%s/ref_00%i.nii.gz',temp_folder,ii-1);
        ref = sprintf('%s/ref_00%i.nii.gz',temp_folder,ii-1);
        Im_gauss = sprintf('%s/im_Gvol%i_00%i.nii.gz',temp_folder,volume,ii-1);
        Im = sprintf('%s/im_vol%i_00%i.nii.gz',temp_folder,volume,ii-1);
        ImPSR = sprintf('%s/PSRslice_00%i.nii.gz',temp_folder,ii-1);
        ImR1 = sprintf('%s/R1slice_00%i.nii.gz',temp_folder,ii-1);
        ImMTR = sprintf('%s/MTRslice_00%i.nii.gz',temp_folder,ii-1);   
    end
    
    if ii < 10
        output_name=sprintf('%s/imv%is0%itoref_2D',temp_folder,volume,ii);
    else
        output_name=sprintf('%s/imv%is%itoref_2D',temp_folder,volume,ii);
    end
    input2D=['antsRegistration --dimensionality 2 --metric CC[' ref_gauss ',' Im_gauss...
        ',1,5,Regular,0.05] --interpolation Linear --transform Rigid[.05] --initial-moving-transform [' ref_gauss...
        ',' Im_gauss ',1] --winsorize-image-intensities [0.01,0.99] --convergence [500x500x1000x1000x1000,1e-6,10] --shrink-factors '...
        '12x8x4x2x1 --smoothing-sigmas 5x4x3x2x1vox --restrict-deformation 1x1x0 -u --output ' output_name];
    unix(input2D);

    if ii < 10
        output_name2=sprintf('%s/imv%is0%itoref_Affine',temp_folder,volume,ii);
    else
        output_name2=sprintf('%s/imv%is%itoref_Affine',temp_folder,volume,ii);
    end
    
    input2DAff=['antsRegistration --dimensionality 2 --metric CC[' ref_gauss ',' Im_gauss...
        ',1,5,Regular,0.05] --interpolation Linear --transform Affine[.05] --initial-moving-transform ' output_name ...
        '0GenericAffine.mat --winsorize-image-intensities [0.01, 0.99] --convergence [500x500x1000x1000x1000,1e-6,10] --shrink-factors '...
        '12x8x4x2x1 --smoothing-sigmas 5x4x3x2x1vox --restrict-deformation 1x1x0 -u --output ' output_name2];
    %produces output_name0GenericAffine.mat
    
    unix(input2DAff);
    
    if ii < 10
        output_nonrigid=sprintf('%s/imv%is0%itoref_NR',temp_folder,volume,ii);
    else
        output_nonrigid=sprintf('%s/imv%is%itoref_NR',temp_folder,volume,ii);
    end
    inputNonRig=['antsRegistration --dimensionality 2 --metric CC[' ref_gauss ',' Im_gauss ...
        ',1,5,Regular,0.05] --interpolation Linear --transform SyN[.05,3,0] --initial-moving-transform ' output_name2 ...
        '0GenericAffine.mat --winsorize-image-intensities [0.01,0.99] --convergence [500x500x1000x1000x1000,1e-6,10] --shrink-factors' ...
        ' 10x6x4x2x1 --smoothing-sigmas 5x3x2x1x0vox -u --restrict-deformation 1x1x0 --output ' output_nonrigid];
    %produces output_name1InverseWarp.nii.gz, output_name1Warp.nii.gz
    unix(inputNonRig);
    
    %Apply to PSR
    if ii < 10
        output_final_path= sprintf('%s/PSRs0%itoref_final.nii.gz',temp_folder,ii);
    else
        output_final_path= sprintf('%s/PSRs%itoref_final.nii.gz',temp_folder,ii);
    end
    inputApplyT=['antsApplyTransforms --dimensionality 2 -i ' ImPSR ' -o ' output_final_path ' -r ' ...
        ref ' -t ' output_nonrigid '1Warp.nii.gz ' output_nonrigid '0GenericAffine.mat -n Linear'];
    unix(inputApplyT);
    
    %Apply to MT
    if ii < 10
        output_final_path= sprintf('%s/ims0%itoref_final.nii.gz',temp_folder,ii);
    else
        output_final_path= sprintf('%s/ims%itoref_final.nii.gz',temp_folder,ii);
    end
    inputApplyT=['antsApplyTransforms --dimensionality 2 -i ' Im ' -o ' output_final_path ' -r ' ...
        ref ' -t ' output_nonrigid '1Warp.nii.gz ' output_nonrigid '0GenericAffine.mat -n Linear'];
    unix(inputApplyT);

    %Apply to R1obs
    if ii < 10
        output_final_path= sprintf('%s/R1s0%itoref_final.nii.gz',temp_folder,ii);
    else
        output_final_path= sprintf('%s/R1s%itoref_final.nii.gz',temp_folder,ii);
    end
    inputApplyT=['antsApplyTransforms --dimensionality 2 -i ' ImR1 ' -o ' output_final_path ' -r ' ...
        ref ' -t ' output_nonrigid '1Warp.nii.gz ' output_nonrigid '0GenericAffine.mat -n Linear'];
    unix(inputApplyT);
    
    %Apply to MTR
    if ii < 10
        output_final_path= sprintf('%s/MTRs0%itoref_final.nii.gz',temp_folder,ii);
    else
        output_final_path= sprintf('%s/MTRs%itoref_final.nii.gz',temp_folder,ii);
    end
    inputApplyT=['antsApplyTransforms --dimensionality 2 -i ' ImMTR ' -o ' output_final_path ' -r ' ...
        ref ' -t ' output_nonrigid '1Warp.nii.gz ' output_nonrigid '0GenericAffine.mat -n Linear'];
    unix(inputApplyT);
end

unix(sprintf('fslmerge -z %s_PSR2ref.nii.gz %s/PSRs*toref_final.nii.gz',...
    ScanNumber,temp_folder));
unix(sprintf('fslmerge -z %s_R1obs2ref.nii.gz %s/R1s*toref_final.nii.gz',...
    ScanNumber,temp_folder));
unix(sprintf('fslmerge -z %s_MT2ref.nii.gz %s/ims*toref_final.nii.gz',...
    ScanNumber,temp_folder));
unix(sprintf('fslmerge -z %s_MTR2ref.nii.gz %s/MTRs*toref_final.nii.gz',...
    ScanNumber,temp_folder));

unix(sprintf('rm -r %s',temp_folder)); 

fprintf(sprintf('Processing for %s complete.', ScanNumber));
%%



