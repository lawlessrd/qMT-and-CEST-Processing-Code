%% CEST Post Processing %% 

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


%% Load in required files

tic;
home=pwd;
out = regexp(home, '\d+', 'match');
ScanNumber = [out{1}];

% Get respiration data from scan phys logs
[actual_resp] = getRespRVT_r55_3TA('n');

imgPath = sprintf('%s/%s_Registration',home,ScanNumber);
cd(imgPath);

e=dir(sprintf('%s_Reference.nii.gz',ScanNumber));
a=dir(sprintf('%s_CEST_Reg.nii.gz',ScanNumber));
b=dir(sprintf('%s_WASSR_Reg.nii.gz',ScanNumber));
c=dir(sprintf('%s_MFA_Reg_to_MFA.nii.gz',ScanNumber));
d=dir(sprintf('%s_B1_Reg_to_MFA.nii.gz',ScanNumber));

Ref_file=e.name;
CEST_file=a.name;
WASSR_file=b.name;
MFAtoMFA_file=c.name;
B1toMFA_file=d.name;

Ref=niftiread(Ref_file);
CEST = niftiread(CEST_file);
WASSR = niftiread(WASSR_file);
B1toMFA=niftiread(B1toMFA_file);
T1toT1=niftiread(MFAtoMFA_file);

cest_new=sprintf('%s_CEST_PostProcessing',ScanNumber);
mkdir(cest_new);


%MFA Parameters
p.T1flip = 30; % MFA flip angle (degrees)
p.MFA = 6; % # of MFA dynamics (FAs)
p.T1TR = 50; % T1 TR (ms)

%%
[APTasym, MTRrex] = CEST_postproc_02102020(CEST,WASSR,actual_resp,cest_new,ScanNumber);

%% Save output as nii
cd(imgPath);
CESTinfo = niftiinfo(sprintf('%s_Prereg/%s_CEST_vol_1.nii.gz',ScanNumber,ScanNumber));
APTvol = sprintf('%s/%s_APTasym2CEST.nii',cest_new,ScanNumber);
niftiwrite(APTasym,APTvol,CESTinfo);
unix(sprintf('gzip %s',APTvol));
MTRvol = sprintf('%s/%s_MTRrex2CEST.nii',cest_new,ScanNumber);
niftiwrite(MTRrex,MTRvol,CESTinfo);
unix(sprintf('gzip %s',MTRvol));

CESTref = sprintf('%s/%s_RefSlice.nii',cest_new,ScanNumber);
%Resample ref to fit CEST
niftiwrite(mean(Ref(:,:,6:9),3),CESTref,CESTinfo);
unix(sprintf('gzip %s',CESTref));

%% Get transform from CEST to mFFE

CESTGvol = sprintf('%s_Prereg/%s_CEST_S0Gvol.nii.gz',ScanNumber,ScanNumber);
CESTvol = sprintf('%s_Prereg/%s_CEST_S0vol.nii.gz',ScanNumber,ScanNumber);
Refvol = CESTref;


% cd(home)
temp_folder = 'ANTs_Reg';
mkdir(temp_folder);

volume = 1;

% Break Reference (mFFE) into individual slices

unix(sprintf('fslsplit %s %s/ref_ -z',...
    Refvol,temp_folder));

%Break outputs into individual slices
unix(sprintf('fslsplit %s %s/APTslice_ -z',...
    APTvol,temp_folder));
unix(sprintf('fslsplit %s %s/MTRslice_ -z',...
    MTRvol,temp_folder));


% Break Image into individual slices
unix(sprintf('fslsplit %s %s/im_vol%i_ -z',...
    CESTvol,temp_folder,volume));

unix(sprintf('fslsplit %s %s/im_Gvol%i_ -z',...
    CESTGvol,temp_folder,volume));

num_slices = length(dir(sprintf('%s/ref_0*.nii.gz',temp_folder)));

for ii=1:num_slices
    fprintf('\tProcessing for slice number: %i\n', ii);
    
    if ii-1 < 10
        ref = sprintf('%s/ref_000%i.nii.gz',temp_folder,ii-1);
        Im_gauss = sprintf('%s/im_Gvol%i_000%i.nii.gz',temp_folder,volume,ii-1);
        Im = sprintf('%s/im_vol%i_000%i.nii.gz',temp_folder,volume,ii-1);
        ImAPT = sprintf('%s/APTslice_000%i.nii.gz',temp_folder,ii-1);
        ImMTR = sprintf('%s/MTRslice_000%i.nii.gz',temp_folder,ii-1);
    else
        ref = sprintf('%s/ref_00%i.nii.gz',temp_folder,ii-1);
        Im_gauss = sprintf('%s/im_Gvol%i_00%i.nii.gz',temp_folder,volume,ii-1);
        Im = sprintf('%s/im_vol%i_00%i.nii.gz',temp_folder,volume,ii-1);
        ImAPT = sprintf('%s/APTslice_00%i.nii.gz',temp_folder,ii-1);
        ImMTR = sprintf('%s/MTRslice_00%i.nii.gz',temp_folder,ii-1);
    end
    
    if ii < 10
        output_name=sprintf('%s/imv%is0%itoref_2D',temp_folder,volume,ii);
    else
        output_name=sprintf('%s/imv%is%itoref_2D',temp_folder,volume,ii);
    end
    input2D=['antsRegistration --dimensionality 2 --metric CC[' ref ',' Im...
        ',1,5,Regular,0.05] --interpolation Linear --transform Rigid[.05] --initial-moving-transform [' ref...
        ',' Im ',1] --winsorize-image-intensities [0.01,0.99] --convergence [500x500x1000x1000x1000,1e-6,10] --shrink-factors '...
        '12x8x4x2x1 --smoothing-sigmas 5x4x3x2x1vox --restrict-deformation 1x1x0 -u --output ' output_name];
    unix(input2D);

    if ii < 10
        output_name2=sprintf('%s/imv%is0%itoref_Affine',temp_folder,volume,ii);
    else
        output_name2=sprintf('%s/imv%is%itoref_Affine',temp_folder,volume,ii);
    end
    
    input2DAff=['antsRegistration --dimensionality 2 --metric CC[' ref ',' Im...
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
    inputNonRig=['antsRegistration --dimensionality 2 --metric CC[' ref ',' Im ...
        ',1,5,Regular,0.05] --interpolation Linear --transform SyN[.05,3,0] --initial-moving-transform ' output_name2 ...
        '0GenericAffine.mat --winsorize-image-intensities [0.01,0.99] --convergence [500x500x1000x1000x1000,1e-6,10] --shrink-factors' ...
        ' 10x6x4x2x1 --smoothing-sigmas 5x3x2x1x0vox -u --restrict-deformation 1x1x0 --output ' output_nonrigid];
    %produces output_name1InverseWarp.nii.gz, output_name1Warp.nii.gz
    unix(inputNonRig);
    
    %Apply to APT
    if ii < 10
        output_final_path= sprintf('%s/APTs0%itoref_final.nii.gz',temp_folder,ii);
    else
        output_final_path= sprintf('%s/APTs%itoref_final.nii.gz',temp_folder,ii);
    end
    inputApplyT=['antsApplyTransforms --dimensionality 2 -i ' ImAPT ' -o ' output_final_path ' -r ' ...
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

movefile('ANTs_Reg/MTRs01toref_final.nii.gz', sprintf('%s/%s_MTRrex2ref.nii.gz',cest_new,ScanNumber));
movefile('ANTs_Reg/APTs01toref_final.nii.gz', sprintf('%s/%s_APT2ref.nii.gz',cest_new,ScanNumber));

%% Calculate T1 maps


for s = 1: 11
    
    [row, col] = find(mask(:,:,s));
    
    fprintf('Slice %g: \n', s)
    fprintf('Number of voxels: %g\n',length(row));
    
    for ii = drange(1:length(row))
        
        if mod(ii/1000,1) == 0
            fprintf('%g/%g: \n',ii, length(row));
        end
        % MFA - R1obs - (Ke's Code)
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

niftiwrite(R1obsPrereg,sprintf('%s_R12R1.nii',ScanNumber));
unix(sprintf('gzip %s_R12R1.nii',ScanNumber));

%% Move R1obs to mFFE space
%% Resample PSR, MT to 14 slices Save results

antsPath = '/Users/lawlesrd/install/bin/';

Refvol = sprintf('%s_Prereg/%s_mFFE_Gvol_1.nii.gz',ScanNumber,ScanNumber);
T1vol = sprintf('%s_R12R1.nii',ScanNumber);

R1r = sprintf('%s_Prereg/%s_R12R1_resample.nii.gz',ScanNumber,ScanNumber);

RefInfo = niftiinfo(Refvol);

% cd(home)
temp_folder = 'ANTs_Reg';
mkdir(temp_folder);

% Resample to match mFFE
% (150x150x14)

unix(['sct_resample -i ' T1vol ' -vox ' num2str(RefInfo.ImageSize(1)) ...
    'x' num2str(RefInfo.ImageSize(2)) 'x' num2str(RefInfo.ImageSize(3)) ...
    ' -x nn -o ' R1r ]);

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
    else
        ref_gauss = sprintf('%s/ref_00%i.nii.gz',temp_folder,ii-1);
        ref = sprintf('%s/ref_00%i.nii.gz',temp_folder,ii-1);
        Im_gauss = sprintf('%s/im_Gvol%i_00%i.nii.gz',temp_folder,volume,ii-1);
        Im = sprintf('%s/im_vol%i_00%i.nii.gz',temp_folder,volume,ii-1);
        ImPSR = sprintf('%s/PSRslice_00%i.nii.gz',temp_folder,ii-1);
        ImR1 = sprintf('%s/R1slice_00%i.nii.gz',temp_folder,ii-1);
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
end

unix(sprintf('fslmerge -z %s_PSR2ref.nii.gz %s/PSRs*toref_final.nii.gz',...
    ScanNumber,temp_folder));
unix(sprintf('fslmerge -z %s_R1obs2ref.nii.gz %s/R1s*toref_final.nii.gz',...
    ScanNumber,temp_folder));
unix(sprintf('fslmerge -z %s_MT2ref.nii.gz %s/ims*toref_final.nii.gz',...
    ScanNumber,temp_folder));

unix(sprintf('rm -r %s',temp_folder)); 

fprintf(sprintf('Processing for %s complete.', ScanNumber));

%% AREX

R1map = niftiread(sprintf('%s_R1obs2ref.nii.gz',ScanNumber));

cd(cest_new)
MTRrexRef = niftiread(sprintf('%s/%s_MTRrex2ref.nii.gz',cest_new,ScanNumber));

AREX = MTRrexRef.*mean(R1map(:,:,6:9),3);

% figure(); subplot(1,2,1); imagesc(APTasym); caxis([0 0.5]); ylabel('APT'); colorbar; axis image;
% subplot(122); imagesc(AREX);caxis([0 1.5]); ylabel('AREX'); colorbar; axis image;

niftiwrite(AREX,sprintf('%s_AREX2ref.nii',ScanNumber),CESTinfo);
unix(sprintf('gzip %s_AREX2ref.nii',ScanNumber));
