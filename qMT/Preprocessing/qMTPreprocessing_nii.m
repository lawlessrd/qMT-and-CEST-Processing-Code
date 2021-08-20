% %% qMT Preprocessing Code
% %Load required PAR/REC files
% 
setenv('FSLOUTPUTTYPE','NIFTI_GZ')

    %%
    scanDir = pwd;
    
    out = regexp(scanDir, '\d+', 'match');
    ScanNumber = out{1};
    ScanName = sprintf('Smith_%s',out{1});
    
    a=dir('Smith*WIP_mFFE_0.65.*.nii');
    b=dir('Smith*PulseMT_16_dyn*.nii');
    c=dir('Smith*MFA_T1*.nii');
    d=dir('*B0*Map*.nii');
    
    e=dir('*B1*Map*_e1.nii');
    f=dir('*B1*Map*_e2.nii');
    
    B1anat = e.name;
    B1mag = f.name;
    B1fix( ScanNumber,B1anat,B1mag )
    B0fix(ScanName,d(1).name,d(2).name);
    
    
    g=dir('Smith_*_B1fix*.nii');
    d = dir('Smith_*_B0fix*.nii');
    
    Ref=a.name;
    MT=b.name;
    MFA=c(end).name;
    B0=d.name;
    B1=g.name;
    %%
    Spine_ANTs_MTReg_VU_NoRef(ScanNumber,Ref,MT,B1,B0,MFA,1)

