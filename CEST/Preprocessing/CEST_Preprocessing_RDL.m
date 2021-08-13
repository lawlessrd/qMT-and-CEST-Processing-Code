%% CEST Preprocessing Code
 
%% Set paths

setenv('FSLOUTPUTTYPE','NIFTI_GZ')
Path = getenv('PATH');
if ~contains(Path,'Users/lawlesrd/Documents/fsl/bin')
    setenv('PATH', [Path ':/Users/lawlesrd/Documents/fsl/bin']);  %Edit rdl 4/19/2018 getenv was just PATH
     % setenv('PATH', [Path ':/usr/local/fsl/bin']);  %home version
end

if ~contains(Path,'/Users/lawlesrd/install/bin/')
    setenv('PATH', [Path ':/Users/lawlesrd/install/bin/']);  %Edit rdl 4/19/2018 getenv was just PATH
      %setenv('PATH', [Path ':/opt/ants/bin/']);  %home version
end
clear PATH

%% Load required files
home = pwd;

out = regexp(home, '\d+', 'match');
ScanNumber = ['Smith_' out{1}];

a=dir('*WIP_mFFE_0.65*.nii');
b=dir('*150ms_5NSA*.nii');
c=dir('*WASSR*.nii');
d=dir('*MFA_T1*.nii');
e=dir('*B1*Map*_e1.nii');
f=dir('*B1*Map*_e2.nii');

B1anat = e.name;
B1mag = f.name;
B1fix( ScanName,B1anat,B1mag )

g=dir('Smith*B1fix*.nii');

Ref=a(end).name;
CEST = b.name;
WASSR = c.name;
MFA=d(end).name;
B1=g.name;

%% Run registration function
Spine_ANTs_CEST(ScanNumber,Ref,B1,MFA,WASSR,CEST,1)
