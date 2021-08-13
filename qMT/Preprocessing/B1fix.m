function [  ] = B1fix( ScanName,B1anat,B1mag )
%Fix gStudy export of B1 images

B1a = vuOpenImage(B1anat);
B1b = vuOpenImage(B1mag);

B1size = size(B1a.Data);

Data = zeros([B1size(1),B1size(2),B1size(3),2]);

Data(:,:,:,1) = imrotate(B1a.Data(:,:,:,1),90);
Data(:,:,:,2) = imrotate(B1b.Data,90);

% B1a = Data;

niftiwrite(Data,sprintf('%s_B1fixed.nii',ScanName))


end

