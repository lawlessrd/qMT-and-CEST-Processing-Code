function Reg_Prep_vu_B1nii(Image, Reference, scanID,ImageType,varargin)
%% Load Image Data
%% Reg Prep: Prepare data for Registration
%
% Reg_Prep(Image, Reference, scanID,ImageType {,ImagePath, [resx, resy],detrend_factor,crop,slice_average})
%
% This function will take an image and resize it/change its FOV to match
% the reference image.  It will then save that image as a .nii.gz image
%
% Image - Image being prepared and saved for Registration
%
% Reference - Reference image that Image is being compared against
%
% scanID - Identifier for scan information to be saved as
%
% ImageType - Label to describe what type of scan is being performed
% IMPORTANT --- If preparing B1 or B0 data, use B1 or B0:
%           B1 - Will pull Magnitude and B1 image
%           B0 - Will pull Magnitude and B0 image
%
%
%
% ImagePath - Optional Path to specify alternate file location
%
% [resx, resy] - Optional Parameters to specify resolution other than that
%                of the Reference
%              - Useful when not wanting to significantly increase
%                resolution to prevent image distortions due to smoothing
%
% detrend_factor - A factor that detrends the data for easier registration
%   Default values that can be used:
%       MFA_detrend = [1.8287    1.6001    1.2828    1.1001    1.0000
%       1.1725]
%
% Author: Alex Smith

% Change log:
% -Created: 20151203
% -Added ImagePath Option: 20160209

% Check if scanID is a character
if ~ischar(scanID)
    scanID = num2str(scanID);
end

% Check if make_nii exists
if isempty(which('make_nii'))
    error('make_nii not in path, add function to path and rerun');
end

if nargin > 4 && ~isempty(varargin{1})
    ImagePath = varargin{1};
else
    ImagePath = pwd;
end

% % Import Image - With conditionals for removing extraneous data
Im = niftiread(sprintf('%s/%s',ImagePath,Image));
ImInfo = niftiinfo(Image);
RefInfo=niftiinfo(Reference);
%Im = load_untouch_nii(sprintf('%s/%s',ImagePath,Image)); %RDL



if strcmpi(ImageType,'B0')

    x = imrotate(squeeze(Im(:,:,:,:)),90);
    x = fliplr(x);
    ImInfo.ImageSize(1) = 128;
    ImInfo.PixelDimensions(1) = 1.25;
    
elseif strcmpi(ImageType,'B1')
    % Magnitude and B1 Image
    %ROTATE ADDED FOR GSTUDY DATA
    x = imrotate(squeeze(Im(:,:,:,:)),90);
    x = fliplr(x);
    %x = squeeze(Im(:,:,:,:));
    %ADDED FOR GSTUDY DATA
    ImInfo.ImageSize(1) = 128;
    ImInfo.PixelDimensions(1) = 1.25;
   
elseif strcmpi(ImageType,'CEST') || strcmpi(ImageType,'WASSR')
    if contains(Image,'.nii') == 1
        Im = imrotate(Im,90);
    end
    
    x=Im;
    %x = reshape(Im, [size(Im,1), size(Im,2),1,size(Im,3)]);
elseif strcmpi(ImageType,'mFFE')
    if contains(Image,'.nii') == 1
        Im = imrotate(Im,90);
    end
    x=squeeze(Im);
else
    % Regular Data
    if contains(Image,'.nii') == 1
        Im = imrotate(Im,90);
    end
    x=squeeze(Im);
end


%% Resize Image to Reference Image

% Load Reference Image and Acquire Dimensions
Ref = niftiread(sprintf('%s/%s',ImagePath,Reference));

rmax = size(Ref,1);
cmax = size(Ref,2);
zmax = size(x,3);
tmax = size(x,4);

if nargin > 5
    if ~isempty(varargin{2})
        f = varargin{2};
        rmax = f(1);
        cmax = f(2);
        clear f
    end
end


fprintf('Resize %s to: [%i,%i]\n',ImageType,rmax,cmax);


% Ref_FOV = str2double(strsplit(Ref.Parms.pardef.FOV_ap_fh_rl_mm,' '));
% Ref_FOV = Ref.Parms.fov(1);
Ref_FOV = RefInfo.ImageSize(1).*RefInfo.PixelDimensions(1); %RDL
% Im_FOV = str2double(strsplit(Im.info.pardef.FOV_ap_fh_rl_mm,' '));
% Im_FOV = Im.Parms.fov(1);

if strcmpi(ImageType,'MT')
    Im_res = 0.625;
elseif strcmpi(ImageType,'MFA')
    Im_res = 0.625;
    
else
    Im_res = ImInfo.PixelDimensions(1);
end

Im_FOV = ImInfo.ImageSize(1).*Im_res;





% Need to include conditional in case FOV is not the same for Ref and Im
if Ref_FOV < Im_FOV
    % This will reduce Im to the FOV of the Ref image
    Im_red = round((Im_FOV - Ref_FOV)/(2*Im_res));
    Im_orig = x(Im_red+1:end-Im_red,Im_red+1:end-Im_red,:,:);
    
    % Need to perform resize so both Ref and Im are at same voxel size
    tmp_Im = zeros(rmax,cmax,zmax,tmax);
    if ndims(x) > 3
        for ii = 1:zmax
            for jj = 1:tmax
                tmp_Im(:,:,ii,jj) = imresize(Im_orig(:,:,ii,jj),[rmax,cmax]);
            end
        end
    else
        for ii = 1:zmax
            tmp_Im(:,:,ii) = imresize(Im_orig(:,:,ii),[rmax,cmax]);
        end
    end
    
    new_Im = tmp_Im;
    
elseif Ref_FOV > Im_FOV
    % This will zero pad Im to the same FOV of Ref
    Im_inc = round((Ref_FOV - Im_FOV)/(2*Im_res));
    Im_orig = padarray(x,[Im_inc Im_inc]);
    
    % Need to perform resize so both Ref and Im are at same voxel size
    tmp_Im = zeros(rmax,cmax,zmax,tmax);
    %     tmp_Im = zeros(size(x));
    if ndims(x) > 3
        for ii = 1:zmax
            for jj = 1:tmax
                tmp_Im(:,:,ii,jj) = imresize(Im_orig(:,:,ii,jj),[rmax,cmax]);
            end
        end
    else
        for ii = 1:zmax
            tmp_Im(:,:,ii) = imresize(Im_orig(:,:,ii),[rmax,cmax]);
        end
    end
    
    new_Im = tmp_Im;
elseif Ref_FOV == Im_FOV
    % This will happen in Ref and Im have same FOV
    
    % Include conditional for 4-D Data
    if ndims(x) > 3
        tmax = size(x,4);
        new_Im = zeros(rmax,cmax,zmax,tmax);
        for ii = 1:size(x,3)
            for jj = 1:size(x,4)
                x_slice=x(:,:,ii,jj);
                x_resized=imresize(x_slice, [rmax,cmax]);
                new_Im(:,:,ii,jj)=x_resized;
            end
        end
    else
        new_Im = zeros(rmax,cmax,zmax);
        for ii = 1:size(x,3)
            x_slice=x(:,:,ii);
            x_resized=imresize(x_slice, [rmax,cmax]);
            new_Im(:,:,ii)=x_resized;
        end
    end
else
    error('FOV''s for Im and Ref were not calculated');
end



%% Crop!

% Set Crop Boundaries
if nargin > 7
    rad = varargin{4};
else
    rad = round(rmax*0.3,-1);
end
Crop1 = rad;

% Set Window Parameters
if isempty(which('Window_Size'))
    Wins = [402 89 800 1000];
else
    Wins = Window_Size;
end

figure(10);
fprintf('Crop %s\n',ImageType);
for ii=1:size(new_Im,3)
    I=new_Im(:,:,ii,1);
    imagesc(I);
    colormap gray
    set(gcf,'Position',Wins);
    axis image;
    grid on;
    set(gca,'yticklabel',[])
    set(gca,'xticklabel',[])
    title(sprintf('Select Center of Spinal Cord for: %s',ImageType));
    [a,b]=ginput(1);
    xmin=round(a-(Crop1/2));
    ymin=round(b-(Crop1/2));
    center_mffe(:,ii)=[a; b];
    
    for volume=1:size(new_Im,4)
        I_ref=new_Im(:,:,ii,volume);
        Icrop_Im(:,:,ii,volume)=imcrop(I_ref,[xmin ymin Crop1-1 Crop1-1]);
    end
end



close(10);  
pause(0.5);

%% Average Over Slices
% Useful if comparing to slices with different volumes, such as for the
% CEST data colleced (where it is one large slice)

if nargin > 8
    av_Icrop = zeros(size(Icrop_Im,1),size(Icrop_Im,2),size(Icrop_Im,4));
    for volume = 1:size(new_Im,4)
        tmp = Icrop_Im(:,:,:,volume);
        av_Icrop(:,:,volume) = mean(tmp(:,:,varargin{5}),3);
    end
    Icrop_Im = av_Icrop;
end


%% Save Images

% Create Directories
% Create Overall Directory
if ~exist(sprintf('%s_Registration',scanID),'dir')
    mkdir(sprintf('%s_Registration',scanID));
end

if ~exist(sprintf('%s_Registration/%s_Prereg',scanID,scanID),'dir')
    mkdir(sprintf('%s_Registration/%s_Prereg',scanID,scanID));
end


% Save each dynamic in the directory as a nii.gz
if nargin > 8
    for tt = 1:tmax
        fprintf('Saving %s Volume %i out\n',ImageType,tt);
        niistruct2=(Icrop_Im(:,:,tt));
        
        filename=sprintf('%s_Registration/%s_Prereg/%s_%s_vol_%i.nii',scanID,scanID,scanID,ImageType,tt);
        niftiwrite(niistruct2, filename) %saves the new cropped image as a .nii.gz file named "filename"
        gzip(filename);
        delete(filename);
    end        
else
    for tt = 1:tmax
        fprintf('Saving %s Volume %i out\n',ImageType,tt);
        niistruct2=(Icrop_Im(:,:,:,tt));
        
        filename=sprintf('%s_Registration/%s_Prereg/%s_%s_vol_%i.nii',scanID,scanID,scanID,ImageType,tt);
        niftiwrite(niistruct2, filename) %saves the new cropped image as a .nii.gz file named "filename"
        gzip(filename);
        delete(filename);
    end
end




% If including a detrending function, detrend the Gaussian Data
if nargin > 6 && ~isempty(varargin{3})
    detrend_factor = varargin{3};
    
    [rcrop_max,ccrop_max,~,~] = size(Icrop_Im);
    for ii = 1:rcrop_max
        for jj = 1:ccrop_max
            for kk = 1:zmax
                Icrop_Im(ii,jj,kk,:) = squeeze(Icrop_Im(ii,jj,kk,:)).*detrend_factor';
            end
        end
    end
end


% Save each Gaussian Transformed dynamic in the directory as a nii.gz
if nargin > 7
    % Create Gaussian kernel (will multiply images by this before registration)
    kernel = fspecial('gaussian',Crop1,Crop1/2);
    for tt = 1:tmax
        fprintf('Saving %s Gauss Volume %i out\n',ImageType,tt);
        niistruct2=(Icrop_Im(:,:,tt).*kernel);
        
        filename=sprintf('%s_Registration/%s_Prereg/%s_%s_Gvol_%i.nii',scanID,scanID,scanID,ImageType,tt);
        niftiwrite(niistruct2, filename) %saves the new cropped image as a .nii.gz file named "filename"
        gzip(filename);
        delete(filename);
    end
else
    % Create Gaussian kernel (will multiply images by this before registration)
    zmax = size(Icrop_Im,3);
    kernel = repmat(fspecial('gaussian',Crop1,Crop1/2),[1 1 zmax]);
    for tt = 1:tmax
        fprintf('Saving %s Gauss Volume %i out\n',ImageType,tt);
        niistruct2=(Icrop_Im(:,:,:,tt).*kernel);
        
        filename=sprintf('%s_Registration/%s_Prereg/%s_%s_Gvol_%i.nii',scanID,scanID,scanID,ImageType,tt);
        niftiwrite(niistruct2, filename) %saves the new cropped image as a .nii.gz file named "filename"
        gzip(filename);
        delete(filename);
    end
end

if length(Image) == length(Reference)
    if sum(Image == Reference) == length(Image)
        fprintf('Saving Reference Full Volume out \n');
        niistruct2 = (Icrop_Im(:,:,:,end));
        
        filename = sprintf('%s_Registration/%s_Reference.nii',scanID,scanID);
        niftiwrite(niistruct2, filename) %saves the new cropped image as a .nii.gz file named "filename"
        gzip(filename);
        delete(filename);
    end
end


%Save center of each image type for segmentation
cd(sprintf('%s_Registration/%s_Prereg',scanID,scanID));
save(sprintf('center_%s.mat',ImageType), 'center_mffe');
cd ..;
cd ..;
end
