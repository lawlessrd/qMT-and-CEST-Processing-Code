function Reg_ANTs_home_NoRef(scanID,ImageType)
%
% Assumes FSL and ANTS are installed on computer.
% For FSL, see: http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/
% For ANTS, see: http://stnava.github.io/ANTs/
% rdl: ANTs changed to ants   04/19/2018

fprintf('\n--------------------------------------- \n');
fprintf('Beginning Registration of %s to MT. \n',ImageType);
fprintf('--------------------------------------- \n');

% Check if loadPARREC exists
if isempty(which('loadPARREC'))
    error('loadPARREC not in path, add function to path and rerun');
end

% Corrects PATH so fsl can be called without listing full path every time
Path = getenv('PATH');
if ~contains(Path,'/Users/lawlesrd/Documents/fsl/bin')
    setenv('PATH', [Path ':/Users/lawlesrd/Documents/fsl/bin']);  %Edit rdl 4/19/2018 getenv was just PATH
      %setenv('PATH', [Path ':/usr/local/fsl/bin']);  %home version
end
clear PATH

Path = getenv('PATH');
if ~contains(Path,'/Users/lawlesrd/install/bin/')
    %setenv('PATH', [Path ':/Users/lawlesrd/Documents/fsl/bin']);  %Edit rdl 4/19/2018 getenv was just PATH
      setenv('PATH', [Path ':/Users/lawlesrd/install/bin/']);  %home version
end
clear PATH

% Check if scanID is a character
if ~ischar(scanID)
    scanID = num2str(scanID);
end

num_dyn = length(dir(sprintf('%s_Registration/%s_Prereg/%s_MT_G*',scanID,scanID,scanID)));

temp_folder='ANTs_Reg';

% Create temp folder for registration
mkdir(temp_folder);

if strcmpi(ImageType,'B1')
    im_reg_B1(scanID,ImageType,temp_folder);
    
%     unix(sprintf('cp %s/imv%itoref.nii.gz %s_%s_Reg.nii.gz',temp_folder,2,scanID,ImageType));

elseif strcmpi(ImageType,'B0')
    im_reg_B0(scanID,ImageType,temp_folder);
    
    unix(sprintf('cp %s/imv%itoref.nii.gz %s_Registration/%s_%s_Reg.nii.gz',temp_folder,2,scanID,scanID,ImageType));
    
elseif strcmpi(ImageType,'MFA') || strcmpi(ImageType,'T115') || strcmpi(ImageType,'T160') || strcmpi(ImageType,'T1')
        im_reg_MFA(scanID,ImageType,num_dyn,temp_folder);
    
elseif strcmpi(ImageType,'CEST') || strcmpi(ImageType,'WASSR')
        im_reg_CEST(scanID,ImageType,num_dyn,temp_folder);
        
else
    if strcmpi(ImageType,'MT') ||  strcmpi(ImageType,'mFFE') || strcmpi(ImageType,'ihMT')
        im_reg_MT(scanID,ImageType,num_dyn,temp_folder);
    
    else
        error('Need to find parameters for new ImageType: %s',ImageType);
    end
    
    fslm_in = [];
    for volume = 1:num_dyn
        fslm_in = sprintf('%s %s',fslm_in,sprintf('%s/imv%itoref.nii.gz',temp_folder,volume));
    end
    
    unix(sprintf('fslmerge -t %s_Registration/%s_%s_Reg %s',scanID,scanID,ImageType,fslm_in));
end

unix(sprintf('rm -r %s',temp_folder));

fprintf('\n------------------------------------------ \n');
fprintf('%s Registration to mFFE Completed. \n',ImageType);
fprintf('------------------------------------------ \n');
end

%% Function for MT weighted Images
function im_reg_MT(scanID,ImageType,num_dyn,temp_folder)

for volume=1:num_dyn
    
    fprintf('\nProcessing for volume number: %i \n', volume);
    
    % Break Reference (mFFE) into individual slices
    if volume == 1
        unix(sprintf('fslsplit %s_Registration/%s_Prereg/%s_MT_vol_1.nii.gz %s/ref_ -z',...
            scanID,scanID,scanID,temp_folder));
        
        unix(sprintf('fslsplit %s_Registration/%s_Prereg/%s_MT_Gvol_1.nii.gz %s/ref_G_ -z',...
            scanID,scanID,scanID,temp_folder));
    end
    
    % Break Image into individual slices
    unix(sprintf('fslsplit %s_Registration/%s_Prereg/%s_%s_vol_%i.nii.gz %s/im_vol%i_ -z',...
        scanID,scanID,scanID,ImageType,volume,temp_folder,volume));
    
    unix(sprintf('fslsplit %s_Registration/%s_Prereg/%s_%s_Gvol_%i.nii.gz %s/im_Gvol%i_ -z',...
        scanID,scanID,scanID,ImageType,volume,temp_folder,volume));
    
    num_slices = length(dir(sprintf('%s/ref_0*.nii.gz',temp_folder)));
    
    for ii=1:num_slices
        fprintf('\tProcessing for slice number: %i\n', ii);
        
        if ii-1 < 10
            ref_gauss = sprintf('%s/ref_G_000%i.nii.gz',temp_folder,ii-1);
            ref = sprintf('%s/ref_000%i.nii.gz',temp_folder,ii-1);
            Im_gauss = sprintf('%s/im_Gvol%i_000%i.nii.gz',temp_folder,volume,ii-1);
            Im = sprintf('%s/im_vol%i_000%i.nii.gz',temp_folder,volume,ii-1);
        else
            ref_gauss = sprintf('%s/ref_00%i.nii.gz',temp_folder,ii-1);
            ref = sprintf('%s/ref_00%i.nii.gz',temp_folder,ii-1);
            Im_gauss = sprintf('%s/im_Gvol%i_00%i.nii.gz',temp_folder,volume,ii-1);
            Im = sprintf('%s/im_vol%i_00%i.nii.gz',temp_folder,volume,ii-1);
        end
        
        if ii < 10
            output_name=sprintf('%s/imv%is0%itoref_2D',temp_folder,volume,ii);
        else
            output_name=sprintf('%s/imv%is%itoref_2D',temp_folder,volume,ii);
        end
        input2D=[' antsRegistration --dimensionality 2 --metric CC[' ref_gauss ',' Im_gauss...
            ',1,5,Regular,0.05] --interpolation Linear --transform Rigid[.05] --initial-moving-transform [' ref_gauss...
            ',' Im_gauss ',1] --winsorize-image-intensities [0.01,0.99] --convergence [500x500x1000x1000x1000,1e-6,10] --shrink-factors '...
            '12x8x4x2x1 --smoothing-sigmas 5x4x3x2x1vox --restrict-deformation 1x1x0 -u --output ' output_name];
        unix(input2D);
        
        if ii < 10
            output_name2=sprintf('%s/imv%is0%itoref_Affine',temp_folder,volume,ii);
        else
            output_name2=sprintf('%s/imv%is%itoref_Affine',temp_folder,volume,ii);
        end
        
        input2DAff=[' antsRegistration --dimensionality 2 --metric CC[' ref_gauss ',' Im_gauss...
            ',1,5,Regular,0.05] --interpolation Linear --transform Affine[.05] --initial-moving-transform ' output_name ...
            '0GenericAffine.mat --winsorize-image-intensities [0.01, 0.99] --convergence [500x500x1000x1000x1000,1e-6,10] --shrink-factors '...
            '12x8x4x2x1 --smoothing-sigmas 5x4x3x2x1vox --restrict-deformation 1x1x0 -u --output ' output_name2];
        %produces output_name0GenericAffine.mat
        
        unix(input2DAff);
        %[statusAff, cmndAff]=unix(input2DAff);
        %outputstring = cmndAff;
        %fprintf(fileid, '%s\n', outputString);
%         
%         if ii < 10
%             output_nonrigid=sprintf('%s/imv%is0%itoref_NR',temp_folder,volume,ii);
%         else
%             output_nonrigid=sprintf('%s/imv%is%itoref_NR',temp_folder,volume,ii);
%         end
%         inputNonRig=[' antsRegistration --dimensionality 2 --metric CC[' ref_gauss ',' Im_gauss ...
%             ',1,5,Regular,0.05] --interpolation Linear --transform SyN[.05,3,0] --initial-moving-transform ' output_name2 ...
%             '0GenericAffine.mat --winsorize-image-intensities [0.01,0.99] --convergence [500x500x1000x1000x1000,1e-6,10] --shrink-factors' ...
%             ' 10x6x4x2x1 --smoothing-sigmas 5x3x2x1x0vox -u --restrict-deformation 1x1x0 --output ' output_nonrigid];
%         %produces output_name1InverseWarp.nii.gz, output_name1Warp.nii.gz
%         unix(inputNonRig);
        
        if ii < 10
            output_final_path= sprintf('%s/imv%is0%itoref_final.nii.gz',temp_folder,volume,ii);
        else
            output_final_path= sprintf('%s/imv%is%itoref_final.nii.gz',temp_folder,volume,ii);
        end
        inputApplyT=[' antsApplyTransforms --dimensionality 2 -i ' Im ' -o ' output_final_path ' -r ' ...
            ref ' -t ' output_name2 '0GenericAffine.mat -n Linear']; %0GenericAffine
        unix(inputApplyT);
        %[statusNR, cmndNR]=unix(inputApplyT);
        %outputstring=cmndNR;
        %fprintf(fileid, '%s\n', outputString);
    end
    
    unix(sprintf('fslmerge -z %s/imv%itoref.nii.gz %s/imv%is*toref_final.nii.gz',...
        temp_folder,volume,temp_folder,volume));
    
end

end


%% Function for MFA Images
function im_reg_MFA(scanID,ImageType,num_dyn,temp_folder)

fprintf('\nProcessing for %s \n', ImageType);

for ii = num_dyn:-1:2
    in_vol = sprintf('%s_Registration/%s_Prereg/%s_%s_Gvol_%i.nii.gz',scanID,scanID,scanID,ImageType,ii);
    trans_mat = sprintf('%s/%s_T1_vol%i_to_%i.xfm',temp_folder,scanID,ii,ii-1);
    target_vol = sprintf('%s_Registration/%s_Prereg/%s_%s_Gvol_%i.nii.gz',scanID,scanID,scanID,ImageType,ii-1);
    
    input_trans = ['flirt -2D -searchry -5 5 -searchrx -5 5 -searchrz 0 0'...
        ' -interp spline -searchcost normcorr -datatype float -in '...
        in_vol ' -ref ' target_vol ' -omat ' trans_mat];
    
    disp(['Find T1 Coreg Trans Mat for: ' in_vol ' to ' target_vol])
    unix(input_trans);
end

% Find Transform to Register T1 to the first volume
in_vol = sprintf('%s_Registration/%s_Prereg/%s_%s_Gvol_1.nii.gz',scanID,scanID,scanID,ImageType);
trans_mat_1 = sprintf('%s/%s_T1_vol1_to_T1.xfm',temp_folder,scanID);
target_vol = sprintf('%s_Registration/%s_Prereg/%s_MFA_Gvol_1.nii.gz',scanID,scanID,scanID);

input_trans = ['flirt -2D '...
    '-searchry -5 5 -searchrx -5 5 -searchrz 0 0 '...
    '-searchcost normcorr -datatype float -in ' in_vol ' -ref '...
    target_vol ' -omat ' trans_mat_1];

disp(['Find T1 Reg Trans Mat for: ' in_vol ' to ' target_vol])
unix(input_trans);



% Concatenate Each Transform Matrix Together
for ii = 1:num_dyn
    mat_name{ii} = sprintf('%s/%s_T1_vol%i_to_%i.xfm',temp_folder,scanID,ii+1,ii);
end

% Find Transform from Dynamic 2 to mFFE
% Repeats added by rdl   04/20/2018
unix(sprintf('convert_xfm -omat %s/%s_T1_vol2_to_T1.xfm -concat %s/%s_T1_vol1_to_T1.xfm %s', temp_folder,scanID,temp_folder,scanID,mat_name{1}));
disp('Find T1 Reg Trans Mat For: T1 2 to first volume')

unix(sprintf('convert_xfm -omat %s/%s_T1_vol3_to_T1.xfm -concat %s/%s_T1_vol2_to_T1.xfm %s', temp_folder,scanID,temp_folder,scanID,mat_name{2}));
disp('Find T1 Reg Trans Mat For: T1 3 to first volume')

unix(sprintf('convert_xfm -omat %s/%s_T1_vol4_to_T1.xfm -concat %s/%s_T1_vol3_to_T1.xfm %s', temp_folder,scanID,temp_folder,scanID,mat_name{3}));
disp('Find T1 Reg Trans Mat For: T1 4 to first volume')

unix(sprintf('convert_xfm -omat %s/%s_T1_vol5_to_T1.xfm -concat %s/%s_T1_vol4_to_T1.xfm %s', temp_folder,scanID,temp_folder,scanID,mat_name{4}));
disp('Find T1 Reg Trans Mat For: T1 5 to first volume')

unix(sprintf('convert_xfm -omat %s/%s_T1_vol6_to_T1.xfm -concat %s/%s_T1_vol5_to_T1.xfm %s', temp_folder,scanID,temp_folder,scanID,mat_name{5}));
disp('Find T1 Reg Trans Mat For: T1 6 to first volume')

for ii = 1:num_dyn-1
    unix(sprintf('convert_xfm -omat %s/T1_reg_mat_%i_to_T1.xfm -concat %s/%s_T1_vol%i_to_T1.xfm %s',temp_folder,ii+1,temp_folder,scanID,ii,mat_name{ii}));
    disp(['Find T1 Reg Trans Mat For: T1 ' num2str(ii+1) ' to MT']);
end

% Perform Registration of T1 data to MT
for ii = 1:num_dyn
    in_vol = sprintf('%s_Registration/%s_Prereg/%s_%s_vol_%i.nii.gz',scanID,scanID,scanID,ImageType,ii);
    trans_mat = sprintf('%s/%s_T1_vol%i_to_T1.xfm',temp_folder,scanID,ii);
    target_vol = sprintf('%s_Registration/%s_Prereg/%s_MFA_vol_1.nii.gz',scanID,scanID,scanID);
    out_vol = sprintf('%s/imv%itoref.nii.gz',temp_folder,ii);
    
    input_trans = sprintf('flirt -ref %s -in %s -2D -applyxfm -init %s -out %s',target_vol,in_vol,trans_mat,out_vol);
    
    disp(['Register T1 Data for: ' in_vol ' to ' target_vol])
    unix(input_trans);
end


% Merge registered MFA dynamics
fslm_in = [];
for volume = 1:num_dyn
    fslm_in = sprintf('%s %s',fslm_in,sprintf('%s/imv%itoref.nii.gz',temp_folder,volume));
end

out_vol = sprintf('%s_Registration/%s_%s_Reg_to_MFA.nii.gz',scanID,scanID,ImageType);

unix(sprintf('fslmerge -t %s %s',out_vol,fslm_in));

%% Find transform from dynamic 1 of MFA to MT

    for ii = num_dyn:-1:2
        in_vol = sprintf('%s_Registration/%s_Prereg/%s_%s_Gvol_%i.nii.gz',scanID,scanID,scanID,ImageType,ii);
        trans_mat = sprintf('%s/%s_T1_vol%i_to_%i.xfm',temp_folder,scanID,ii,ii-1);
        target_vol = sprintf('%s_Registration/%s_Prereg/%s_%s_Gvol_%i.nii.gz',scanID,scanID,scanID,ImageType,ii-1);
        
        input_trans = ['flirt -2D -searchry -5 5 -searchrx -5 5 -searchrz 0 0'...
            ' -interp spline -searchcost normcorr -datatype float -in '...
            in_vol ' -ref ' target_vol ' -omat ' trans_mat];
        
        disp(['Find T1 Coreg Trans Mat for: ' in_vol ' to ' target_vol])
        unix(input_trans);
    end
    
    % Find Transform to Register T1 to MT
    in_vol = sprintf('%s_Registration/%s_Prereg/%s_%s_Gvol_1.nii.gz',scanID,scanID,scanID,ImageType);
    trans_mat_1 = sprintf('%s_Registration/%s_Prereg/%s_T1_vol1_to_MT.xfm',scanID,scanID,scanID);
    target_vol = sprintf('%s_Registration/%s_Prereg/%s_MT_Gvol_1.nii.gz',scanID,scanID,scanID);
    
    input_trans = ['flirt -2D '...
        '-searchry -5 5 -searchrx -5 5 -searchrz 0 0 '...
        '-searchcost normcorr -datatype float -in ' in_vol ' -ref '...
        target_vol ' -omat ' trans_mat_1];
    
    disp(['Find T1 Reg Trans Mat for: ' in_vol ' to ' target_vol])
    unix(input_trans);
    
   
    
    % Concatenate Each Transform Matrix Together
    for ii = 1:num_dyn
        mat_name{ii} = sprintf('%s/%s_T1_vol%i_to_%i.xfm',temp_folder,scanID,ii+1,ii);
    end
    
    % Find Transform from Dynamic 2 to mFFE
    % Repeats added by rdl   04/20/2018
    unix(sprintf('convert_xfm -omat %s/%s_T1_vol2_to_MT.xfm -concat %s %s', temp_folder,scanID,trans_mat_1,mat_name{1}));
    disp('Find T1 Reg Trans Mat For: T1 2 to MT')
    
    unix(sprintf('convert_xfm -omat %s/%s_T1_vol3_to_MT.xfm -concat %s/%s_T1_vol2_to_MT.xfm %s', temp_folder,scanID,temp_folder,scanID,mat_name{2}));
    disp('Find T1 Reg Trans Mat For: T1 3 to MT')
    
    unix(sprintf('convert_xfm -omat %s/%s_T1_vol4_to_MT.xfm -concat %s/%s_T1_vol3_to_MT.xfm %s', temp_folder,scanID,temp_folder,scanID,mat_name{3}));
    disp('Find T1 Reg Trans Mat For: T1 4 to MT')
    
    unix(sprintf('convert_xfm -omat %s/%s_T1_vol5_to_MT.xfm -concat %s/%s_T1_vol4_to_MT.xfm %s', temp_folder,scanID,temp_folder,scanID,mat_name{4}));
    disp('Find T1 Reg Trans Mat For: T1 5 to MT')    
    
    unix(sprintf('convert_xfm -omat %s/%s_T1_vol6_to_MT.xfm -concat %s/%s_T1_vol5_to_MT.xfm %s', temp_folder,scanID,temp_folder,scanID,mat_name{5}));
    disp('Find T1 Reg Trans Mat For: T1 6 to MT')    
    
    for ii = 1:num_dyn-1
        if  ii == 1
            unix(sprintf('convert_xfm -omat %s/T1_reg_mat_%i_to_MT.xfm -concat %s %s',temp_folder,ii+1,trans_mat_1, mat_name{ii}));
            disp(['Find T1 Reg Trans Mat For: T1 ' num2str(ii+1) ' to MT']);
        else
            unix(sprintf('convert_xfm -omat %s/T1_reg_mat_%i_to_MT.xfm -concat %s/%s_T1_vol%i_to_MT.xfm %s',temp_folder,ii+1,temp_folder,scanID,ii,mat_name{ii}));
            disp(['Find T1 Reg Trans Mat For: T1 ' num2str(ii+1) ' to MT']);
        end
    end
    
    % Perform Registration of T1 data to MT
    for ii = 1:num_dyn
        in_vol = sprintf('%s_Registration/%s_Prereg/%s_%s_vol_%i.nii.gz',scanID,scanID,scanID,ImageType,ii);
        if ii == 1
            trans_mat = trans_mat_1;
        else
            trans_mat = sprintf('%s/%s_T1_vol%i_to_MT.xfm',temp_folder,scanID,ii);
        end
        target_vol = sprintf('%s_Registration/%s_Prereg/%s_MT_vol_1.nii.gz',scanID,scanID,scanID);
        out_vol = sprintf('%s/imv%itoMT.nii.gz',temp_folder,ii);
        
        input_trans = sprintf('flirt -ref %s -in %s -2D -applyxfm -init %s -out %s',target_vol,in_vol,trans_mat,out_vol);
        
        disp(['Register T1 Data for: ' in_vol ' to ' target_vol])
        unix(input_trans);
    end
    
    % Merge registered MFA dynamics
    fslm_in = [];
    for volume = 1:num_dyn
        fslm_in = sprintf('%s %s',fslm_in,sprintf('%s/imv%itoMT.nii.gz',temp_folder,volume));
    end
    
    out_vol = sprintf('%s_Registration/%s_%s_Reg_to_MT.nii.gz',scanID,scanID,ImageType);
    
    unix(sprintf('fslmerge -t %s %s',out_vol,fslm_in));

end

%% Function for B0 Images
function im_reg_B0(scanID,ImageType,temp_folder)

fprintf('\nProcessing for %s \n', ImageType);

% Break Reference (MT) into individual slices
unix(sprintf('fslsplit %s_Registration/%s_Prereg/%s_MT_vol_1.nii.gz %s/ref_ -z',...
    scanID,scanID,scanID,temp_folder));

unix(sprintf('fslsplit %s_Registration/%s_Prereg/%s_MT_Gvol_1.nii.gz %s/ref_G_ -z',...
    scanID,scanID,scanID,temp_folder));

% Break Image into individual slices
unix(sprintf('fslsplit %s_Registration/%s_Prereg/%s_%s_vol_1.nii.gz %s/im_vol1_ -z',...
    scanID,scanID,scanID,ImageType,temp_folder));

unix(sprintf('fslsplit %s_Registration/%s_Prereg/%s_%s_vol_2.nii.gz %s/im_vol2_ -z',...
    scanID,scanID,scanID,ImageType,temp_folder));

unix(sprintf('fslsplit %s_Registration/%s_Prereg/%s_%s_Gvol_1.nii.gz %s/im_Gvol1_ -z',...
    scanID,scanID,scanID,ImageType,temp_folder));

unix(sprintf('fslsplit %s_Registration/%s_Prereg/%s_%s_Gvol_2.nii.gz %s/im_Gvol2_ -z',...
    scanID,scanID,scanID,ImageType,temp_folder));

num_slices = length(dir(sprintf('%s/ref_0*.nii.gz',temp_folder)));

for ii=1:num_slices
    fprintf('\tProcessing for slice number: %i\n', ii);
    
    if ii-1 < 10
        ref_gauss = sprintf('%s/ref_G_000%i.nii.gz',temp_folder,ii-1);
        ref = sprintf('%s/ref_000%i.nii.gz',temp_folder,ii-1);
        Im_gauss = sprintf('%s/im_Gvol1_000%i.nii.gz',temp_folder,ii-1);
        Im = sprintf('%s/im_vol2_000%i.nii.gz',temp_folder,ii-1);
    else
        ref_gauss = sprintf('%s/ref_00%i.nii.gz',temp_folder,ii-1);
        ref = sprintf('%s/ref_00%i.nii.gz',temp_folder,ii-1);
        Im_gauss = sprintf('%s/im_Gvol1_00%i.nii.gz',temp_folder,ii-1);
        Im = sprintf('%s/im_vol2_00%i.nii.gz',temp_folder,ii-1);
    end
    
    if ii < 10
        output_name=sprintf('%s/imv%is0%itoref_2D',temp_folder,1,ii);
    else
        output_name=sprintf('%s/imv%is%itoref_2D',temp_folder,1,ii);
    end
    input2D=[' antsRegistration --dimensionality 2 --metric CC[' ref_gauss ',' Im_gauss...
        ',1,5,Regular,0.05] --interpolation NearestNeighbor --transform Rigid[.05] --initial-moving-transform [' ref_gauss...
        ',' Im_gauss ',1] --winsorize-image-intensities [0.01,0.99] --convergence [500x500x1000x1000x1000,1e-6,10] --shrink-factors '...
        '12x8x4x2x1 --smoothing-sigmas 5x4x3x2x1vox --restrict-deformation 1x1x0 -u --output ' output_name];
    unix(input2D);
    
    if ii < 10
        output_final_path= sprintf('%s/imv%is0%itoref_final.nii.gz',temp_folder,2,ii);
    else
        output_final_path= sprintf('%s/imv%is%itoref_final.nii.gz',temp_folder,2,ii);
    end
    inputApplyT=[' antsApplyTransforms --dimensionality 2 -i ' Im ' -o ' output_final_path ' -r ' ...
        ref ' -t ' output_name '0GenericAffine.mat -n Linear'];
    unix(inputApplyT);
    
end

unix(sprintf('fslmerge -z %s/imv%itoref.nii.gz %s/imv%is*toref_final.nii.gz',...
    temp_folder,2,temp_folder,2));

end

%% Function for B1 Images
function im_reg_B1(scanID,ImageType,temp_folder)

% Using FSL's FLIRT function for the B1 map as ANTS is not well-determined
% for the B1 map.

fprintf('\nProcessing for %s \n', ImageType);

% Find Transform to Register B1 to the MT
in_vol = sprintf('%s_Registration/%s_Prereg/%s_B1_Gvol_1.nii.gz',scanID,scanID,scanID);
trans_mat = sprintf('%s/%s_B1_to_MT.mat',temp_folder,scanID);
target_vol = sprintf('%s_Registration/%s_Prereg/%s_MT_Gvol_1.nii.gz',scanID,scanID,scanID);

input_trans = sprintf(['flirt -2D -searchry 0 0 -searchrx 0 0 '...
    '-searchrz -10 10 -interp spline -searchcost normcorr -datatype float'...
    ' -in %s -ref %s -omat %s'],in_vol,target_vol,trans_mat);

disp(['Find B1 Reg Trans Mat for: ' in_vol ' to ' target_vol])
unix(input_trans);

% Apply the Tranform Matrix to the B1 Map
in_vol = sprintf('%s_Registration/%s_Prereg/%s_B1_vol_2.nii.gz',scanID,scanID,scanID);
trans_mat = sprintf('%s/%s_B1_to_MT.mat',temp_folder,scanID);
target_vol = sprintf('%s_Registration/%s_Prereg/%s_MT_vol_1.nii.gz',scanID,scanID,scanID);
out_vol = sprintf('%s_Registration/%s_B1_Reg_to_MT.nii.gz',scanID,scanID);

input_trans = sprintf('flirt -ref %s -in %s -2D -applyxfm -init %s -out %s',...
    target_vol,in_vol,trans_mat,out_vol);

disp(['Register B1 Data for: ' in_vol ' to ' target_vol])
unix(input_trans);



% Find Transform to Register B1 to the MFA
in_vol = sprintf('%s_Registration/%s_Prereg/%s_B1_Gvol_1.nii.gz',scanID,scanID,scanID);
trans_mat = sprintf('%s/%s_B1_to_MFA.mat',temp_folder,scanID);
target_vol = sprintf('%s_Registration/%s_Prereg/%s_MFA_Gvol_1.nii.gz',scanID,scanID,scanID);

input_trans = sprintf(['flirt -2D -searchry 0 0 -searchrx 0 0 '...
    '-searchrz -10 10 -interp spline -searchcost normcorr -datatype float'...
    ' -in %s -ref %s -omat %s'],in_vol,target_vol,trans_mat);

disp(['Find B1 Reg Trans Mat for: ' in_vol ' to ' target_vol])
unix(input_trans);

% Apply the Tranform Matrix to the B1 Map
in_vol = sprintf('%s_Registration/%s_Prereg/%s_B1_vol_2.nii.gz',scanID,scanID,scanID);
trans_mat = sprintf('%s/%s_B1_to_MFA.mat',temp_folder,scanID);
target_vol = sprintf('%s_Registration/%s_Prereg/%s_MFA_Gvol_1.nii.gz',scanID,scanID,scanID);
out_vol = sprintf('%s_Registration/%s_B1_Reg_to_MFA.nii.gz',scanID,scanID);

input_trans = sprintf('flirt -ref %s -in %s -2D -applyxfm -init %s -out %s',...
    target_vol,in_vol,trans_mat,out_vol);

disp(['Register B1 Data for: ' in_vol ' to ' target_vol])
unix(input_trans);

end
%% Function for MT weighted Images
function im_reg_CEST(scanID,ImageType,num_dyn,temp_folder)

    % Save niftis of MFA, mFFE mean of slices 6:9
    
    
    % Get transform from MFA to CEST (Postprocessing? Wait until they are
    % both in mFFE sapce?
    
    % Save S0 mean as target
    
    
    %% Load in offset list to determine S0 images and make mean S0 for target
    
    if strcmpi(ImageType,'CEST')
        if num_dyn == 37
            rf_A=load('/Users/lawlesrd/Desktop/Research/Research Code/CEST/Post Processing/multipleS0_reduced.txt');
        else
            rf_A=load('/Users/lawlesrd/Desktop/Research/Research Code/CEST/Post Processing/multipleS0.txt');
        end
    else
        if num_dyn == 15
            rf_A=load('/Users/lawlesrd/Desktop/Research/Research Code/CEST/Post Processing/Shortened_WASSR.txt');
        else
            rf_A=load('/Users/lawlesrd/Desktop/Research/Research Code/CEST/Post Processing/WASSR.txt');
        end
    end
    
    S0ind = find(rf_A == 100000);
    
    for vol = 1:length(S0ind)
    
        S0vol(:,:,:,vol) = niftiread(sprintf('%s_%s_Gvol_%i.nii.gz',scanID,ImageType,S0ind(i)));
        
        
    end 
    
    
for volume=1:num_dyn
    
    fprintf('\nProcessing for volume number: %i \n', volume);
    
    % Break Reference (mFFE) into individual slices
    if volume == 1
        unix(sprintf('fslsplit %s_Registration/%s_Prereg/%s_MT_vol_1.nii.gz %s/ref_ -z',...
            scanID,scanID,scanID,temp_folder));
        
        unix(sprintf('fslsplit %s_Registration/%s_Prereg/%s_MT_Gvol_1.nii.gz %s/ref_G_ -z',...
            scanID,scanID,scanID,temp_folder));
    end
    
    % Break Image into individual slices
    unix(sprintf('fslsplit %s_Registration/%s_Prereg/%s_%s_vol_%i.nii.gz %s/im_vol%i_ -z',...
        scanID,scanID,scanID,ImageType,volume,temp_folder,volume));
    
    unix(sprintf('fslsplit %s_Registration/%s_Prereg/%s_%s_Gvol_%i.nii.gz %s/im_Gvol%i_ -z',...
        scanID,scanID,scanID,ImageType,volume,temp_folder,volume));
    
    num_slices = length(dir(sprintf('%s/ref_0*.nii.gz',temp_folder)));
    
    for ii=1:num_slices
        fprintf('\tProcessing for slice number: %i\n', ii);
        
        if ii-1 < 10
            ref_gauss = sprintf('%s/ref_G_000%i.nii.gz',temp_folder,ii-1);
            ref = sprintf('%s/ref_000%i.nii.gz',temp_folder,ii-1);
            Im_gauss = sprintf('%s/im_Gvol%i_000%i.nii.gz',temp_folder,volume,ii-1);
            Im = sprintf('%s/im_vol%i_000%i.nii.gz',temp_folder,volume,ii-1);
        else
            ref_gauss = sprintf('%s/ref_00%i.nii.gz',temp_folder,ii-1);
            ref = sprintf('%s/ref_00%i.nii.gz',temp_folder,ii-1);
            Im_gauss = sprintf('%s/im_Gvol%i_00%i.nii.gz',temp_folder,volume,ii-1);
            Im = sprintf('%s/im_vol%i_00%i.nii.gz',temp_folder,volume,ii-1);
        end
        
        if ii < 10
            output_name=sprintf('%s/imv%is0%itoref_2D',temp_folder,volume,ii);
        else
            output_name=sprintf('%s/imv%is%itoref_2D',temp_folder,volume,ii);
        end
        input2D=[' antsRegistration --dimensionality 2 --metric CC[' ref_gauss ',' Im_gauss...
            ',1,5,Regular,0.05] --interpolation Linear --transform Rigid[.05] --initial-moving-transform [' ref_gauss...
            ',' Im_gauss ',1] --winsorize-image-intensities [0.01,0.99] --convergence [500x500x1000x1000x1000,1e-6,10] --shrink-factors '...
            '12x8x4x2x1 --smoothing-sigmas 5x4x3x2x1vox --restrict-deformation 1x1x0 -u --output ' output_name];
        unix(input2D);
        
        if ii < 10
            output_name2=sprintf('%s/imv%is0%itoref_Affine',temp_folder,volume,ii);
        else
            output_name2=sprintf('%s/imv%is%itoref_Affine',temp_folder,volume,ii);
        end
        
        input2DAff=[' antsRegistration --dimensionality 2 --metric CC[' ref_gauss ',' Im_gauss...
            ',1,5,Regular,0.05] --interpolation Linear --transform Affine[.05] --initial-moving-transform ' output_name ...
            '0GenericAffine.mat --winsorize-image-intensities [0.01, 0.99] --convergence [500x500x1000x1000x1000,1e-6,10] --shrink-factors '...
            '12x8x4x2x1 --smoothing-sigmas 5x4x3x2x1vox --restrict-deformation 1x1x0 -u --output ' output_name2];
        %produces output_name0GenericAffine.mat
        
        unix(input2DAff);
        %[statusAff, cmndAff]=unix(input2DAff);
        %outputstring = cmndAff;
        %fprintf(fileid, '%s\n', outputString);
%         
%         if ii < 10
%             output_nonrigid=sprintf('%s/imv%is0%itoref_NR',temp_folder,volume,ii);
%         else
%             output_nonrigid=sprintf('%s/imv%is%itoref_NR',temp_folder,volume,ii);
%         end
%         inputNonRig=[' antsRegistration --dimensionality 2 --metric CC[' ref_gauss ',' Im_gauss ...
%             ',1,5,Regular,0.05] --interpolation Linear --transform SyN[.05,3,0] --initial-moving-transform ' output_name2 ...
%             '0GenericAffine.mat --winsorize-image-intensities [0.01,0.99] --convergence [500x500x1000x1000x1000,1e-6,10] --shrink-factors' ...
%             ' 10x6x4x2x1 --smoothing-sigmas 5x3x2x1x0vox -u --restrict-deformation 1x1x0 --output ' output_nonrigid];
%         %produces output_name1InverseWarp.nii.gz, output_name1Warp.nii.gz
%         unix(inputNonRig);
        
        if ii < 10
            output_final_path= sprintf('%s/imv%is0%itoref_final.nii.gz',temp_folder,volume,ii);
        else
            output_final_path= sprintf('%s/imv%is%itoref_final.nii.gz',temp_folder,volume,ii);
        end
        inputApplyT=[' antsApplyTransforms --dimensionality 2 -i ' Im ' -o ' output_final_path ' -r ' ...
            ref ' -t ' output_name2 '0GenericAffine.mat -n Linear']; %0GenericAffine
        unix(inputApplyT);
        %[statusNR, cmndNR]=unix(inputApplyT);
        %outputstring=cmndNR;
        %fprintf(fileid, '%s\n', outputString);
    end
    
    unix(sprintf('fslmerge -z %s/imv%itoref.nii.gz %s/imv%is*toref_final.nii.gz',...
        temp_folder,volume,temp_folder,volume));
    
end

end

