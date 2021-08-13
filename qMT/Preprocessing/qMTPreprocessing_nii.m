% %% qMT Preprocessing Code
% %Load required PAR/REC files
% 
setenv('FSLOUTPUTTYPE','NIFTI_GZ')

% Loop through full directory of subjects
home = pwd;
dirc=dir(home);
dirc = dirc(find(cellfun(@isdir,{dirc(:).name})));
dirc(1:2) = [];

for subj = 1:length(dirc)
    
    PatientID = dirc(subj).name;
    PatDir = sprintf('%s/%s', dirc(subj).folder, PatientID);
    cd(PatDir)
    %%
    scanDir = pwd;
    
    out = regexp(scanDir, '\d+', 'match');
    ScanNumber = out{1};
    ScanName = sprintf('Smith_%s',out{1});
    
    a=dir('Smith*WIP_mFFE_SENSE_0.6*.nii');
    b=dir('Smith*PulseMT_16_dyn*.nii');
    c=dir('Smith*MFA_T1*.nii');
    d=dir('B0*Map*.nii.gz');
    
    e=dir('B1*Map*_e1.nii.gz');
    f=dir('B1*Map*_e2.nii.gz');
    
    B0fix(ScanName,d(1).name,d(2).name);
    B1fix(ScanName,e.name,f.name)
    
    g=dir('Smith_*_B1fix*.nii');
    d = dir('Smith_*_B0fix*.nii');
    
    Ref=a.name;
    MT=b.name;
    MFA=c(end).name;
    B0=d.name;
    B1=g.name;
   
    Spine_ANTs_MTReg_VU_NoRef(ScanNumber,Ref,MT,B1,B0,MFA,1)
end
