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


%% Load data
tic;
home=pwd;
out = regexp(home, '\d+', 'match');
ScanNumber = [out{1}];
imgPath = sprintf('%s/%s_Registration',home,ScanNumber);
cd(imgPath);

a=dir(sprintf('%s_MT_Reg_SP.nii.gz',ScanNumber));
b=dir(sprintf('%s_MFA_Reg.nii.gz',ScanNumber));
c=dir(sprintf('%s_B0_Reg.nii.gz',ScanNumber));
d=dir(sprintf('%s_B1_Reg.nii.gz',ScanNumber));
e=dir(sprintf('%s_Reference.nii.gz',ScanNumber));
MT_file=a.name;
MFA_file=b.name;
B0_file=c.name;
B1_file=d.name;
Ref_file=e.name;

qMT_s=load_nii(MT_file);
qMT = qMT_s.img;

B0_s=load_nii(B0_file);
B0 = B0_s.img;

B1_s=load_nii(B1_file);
B1 = B1_s.img;

T1mfa_s=load_nii(MFA_file);
T1mfa = T1mfa_s.img;

Ref_s=load_nii(Ref_file);
Ref = Ref_s.img;

%% Input parameters
% p.TR = 50;
% p.qMTflip = 6;
% p.MT_flip = [820];
% p.pwMT = [20e-3];
% p.deltaMT = [2500 100000];
% p.T1flip = 30;
% p.MFA = 6;
% p.T1TR = 20;

p.TR = 57;
p.qMTflip = 6;
p.MT_flip = [820];
p.pwMT = [20e-3];
p.deltaMT = [2500 100000];
p.T1flip = 30;
p.MFA = 6;
p.T1TR = 20;

%% Normalization

M1 = qMT(:,:,:,:);
% M2 = qMT(:,:,:,9:end);

for jj = 1 : size(M1,3)
    for ii = 1 : size(M1,4)
        M1(:,:,jj,ii) = M1(:,:,jj,ii)./M1(:,:,jj,end);
%         M2(:,:,jj,ii) = M2(:,:,jj,ii)./M2(:,:,jj,end);
    end
end

M1(isinf(M1)) = 0 ;
M1(isnan(M1)) = 0 ;

% M2(isinf(M2)) = 0 ;
% M2(isnan(M2)) = 0 ;

%% FITTING
R1obs = zeros(size(B0));
PSR = zeros(size(B0));
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

% Set mask size
maskSize = 50;

unix(sprintf('sct_create_mask -i %s_MT_Reg.nii.gz -p center -size %0.0f -f cylinder -o %s_mask_%0.0f.nii.gz'...
    ,ScanNumber,maskSize,ScanNumber,maskSize));

mask = niftiread(sprintf('%s_mask_%0.0f.nii.gz',ScanNumber,maskSize));

%%
for s = 1:size(M1,3)
    
    [row, col] = find(mask(:,:,s));
    
    fprintf('Slice %g: \n', s)
    fprintf('Number of voxels: %g\n',length(row));
    
    for ii = drange(1:length(row))
        
        if mod(ii/1000,1) == 0
            fprintf('%g/%g: \n',ii, length(row));
        end
        
        Mn1 = squeeze(M1(row(ii), col(ii), s, : ))';
%         Mn2 = squeeze(M2(row(ii), col(ii), s, : ))';
        Mn(1,:) = Mn1;
%         Mn(2,:) = Mn2;
        
        p.B1 = squeeze(B1(row(ii), col(ii), 23));
        p.Ernst = squeeze(T1mfa(row(ii), col(ii),23,:))';
        p.B0 = squeeze(B0(row(ii), col(ii), 23));
        p.M = Mn;
        [PSR(row(ii), col(ii), s),R1obs(row(ii), col(ii), s),chi2(row(ii), col(ii), s),chi2p(row(ii), col(ii), s),res,resn(row(ii), col(ii), s)] = Analysis_Yarnykh_1pt(p);
%         fprintf('PSR=%g, kba=%g, T2a=%g, T2b=%g, T1obs=%g, T1a=%g, chi2=%g, chi2p=%g, resn=%g, sigma=%g \n',PSR(row(ii), col(ii), s),kba(row(ii), col(ii), s),T2a(row(ii), col(ii), s),T2b(row(ii), col(ii), s),T1obs(row(ii), col(ii), s),T1a(row(ii), col(ii), s),chi2(row(ii), col(ii), s),chi2p(row(ii), col(ii), s),resn(row(ii), col(ii), s),sigma(row(ii), col(ii), s));
    end
    
    MPF = PSR ./ (1 + PSR);
%     kab = PSR .* kba;
    save('qMT_results_SinglePoint');
end

%% Save results

T1_s = qMT_s;
T1_s.img = R1obs;
save_nii(T1_s,strcat('T1obs.nii'));

PSR_s = qMT_s;
PSR_s.img = PSR;
save_nii(PSR_s,strcat('PSR.nii'));

