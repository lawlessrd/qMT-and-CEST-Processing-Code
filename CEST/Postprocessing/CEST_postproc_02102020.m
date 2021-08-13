
function [APTasym, MTRrex] = CEST_postproc_02102020(CEST,WASSR,actual_resp,cest_new,ScanNumber);

maskDir = dir(sprintf('%s_mask_*.nii.gz',ScanNumber));

mask = niftiread(maskDir(1).name);

maskCEST = logical(mask(:,:,7));

cd(cest_new);
A_uncorr=double(CEST);

%% Correcting Resp

A=A_uncorr;

if size(A,4) == 37
    rf_A=load('/Users/lawlesrd/Desktop/Research/Research Code/CEST/Post Processing/multipleS0_reduced.txt');
else
    rf_A=load('/Users/lawlesrd/Desktop/Research/Research Code/CEST/Post Processing/multipleS0.txt');
end

jj=1;
for ii=1:4:size(A,4);
    S0(:,:,:,jj)=A(:,:,:,ii);
    jj=jj+1;
end

%if size(S0) == 37


disp('Interpolating S0')
for slice = 1:size(S0,3);
    disp(strcat('Processing Slice Number:  ',num2str(slice)))
    
    for ii = 1:size(S0,1);
        
        if mod(ii, 20) == 0
            disp(strcat('Processing Row Number:  ',num2str(ii)))
        end
        
        for jj = 1:size(S0,2);
            
            if maskCEST(ii,jj, slice) == 1
                
                dd = double(squeeze(S0(ii,jj,slice,:)));
                xx=1:4:size(A,4);
                new_x=1:size(A,4);
                yy=spline(xx, dd, new_x);
                %%yy=interp1(xx, dd, new_x, 'linear');
                %yy=medfilt1
                S0_fit(ii,jj,slice,:) = yy;
                
            else
                S0_fit(ii,jj,slice,:) = zeros(1,size(A,4));
            end
        end
    end
end

%%
figure,plot(new_x,squeeze(S0_fit(75,75,1,:)),new_x,squeeze(A(75,75,1,:)),'-o'),title(sprintf('Voxel in WM for %s',ScanNumber))


CEST=A./S0_fit;

%% Respiration Correction

for slice = 1:size(A,3);
    disp(strcat('Processing Slice Number:  ',num2str(slice)))
    for ii = 1:size(A,1);
        
        if mod(ii, 20) == 0
            disp(strcat('Processing Row Number:  ',num2str(ii)))
        end
        
        for jj = 1:size(A,2);
            
            if maskCEST(ii,jj, slice) == 1
                resp_norm=actual_resp; %([1:9 26:49])?
                %resp_norm=(actual_resp-min(actual_resp))/range(actual_resp); %scales it from 0 to 1
                bb_norm=squeeze(S0_fit(ii,jj,slice,:))./squeeze(mean(S0_fit(ii,jj,slice,:),4));
                
                CEST_filtered(ii,jj,slice,:)=squeeze(CEST(ii,jj,slice,:)) - resp_norm'*(max(bb_norm)-min(bb_norm));
                clear resp_norm bb_norm
            else
                CEST_filtered(ii,jj,slice,:)=zeros(1,size(A,4));
            end
        end
    end
end

jj=1;
for ii=1:size(rf_A)
    a = 1:4:size(A,4);
    if ismember(ii,a)==0 %remove all of the interspersed S0, checks to see if index is in A, returns only false
        new_rf(jj)=rf_A(ii);
        CESTR(:,:,:,jj)=CEST_filtered(:,:,:,ii);
        jj=jj+1;
    else
        continue
    end
end


clear rf_A;
rf_A=new_rf'.*128;

ref=S0_fit; %%
c = clock;
yy = num2str(c(1)); mm = num2str(c(2)); dd = num2str(c(3)); hh = num2str(c(4)); mmm = num2str(c(5));
save_time = strcat(yy,mm,dd,'_',hh,mmm);


savefile = strcat('Raw_CEST_Data_',save_time,'.mat');
save(savefile, 'A', 'rf_A', 'CESTR', 'S0', 'S0_fit', 'CEST', 'actual_resp', 'maskCEST', 'CEST_filtered')


%% WASSR Calculation %%


    disp('Processing WASSR shift map');
    
    A_wassr=double(WASSR);
    
    if size(A_wassr,4) == 15
        rf_wassr=load('/Users/lawlesrd/Desktop/Research/Research Code/CEST/Post Processing/Shortened_WASSR.txt');
    else
        rf_wassr=load('/Users/lawlesrd/Desktop/Research/Research Code/CEST/Post Processing/WASSR.txt');
    end
    
    %     %Remove intersperse
    jj=1;
    
    if size(A_wassr,4) == 15
        for ii=[1 5 11 15]
            S0_wassr(:,:,:,jj)=A_wassr(:,:,:,ii);
            jj=jj+1;
        end
    else
        for ii=1:4:size(A_wassr,4);
            S0_wassr(:,:,:,jj)=A_wassr(:,:,:,ii);
            jj=jj+1;
        end
    end
    
    ref_wassr=(1/2).*(S0_wassr(:,:,:,1) + S0_wassr(:,:,:,end));
    
    
    for ii = 1:size(A_wassr,4)
        CEST_wassr(:,:,:,ii) = A_wassr(:,:,:,ii)./ref_wassr;
    end
    
    
    jj=1;
    for ii=1:size(rf_wassr)
        if size(A_wassr,4) == 15
            a=[1 5 11 15];
        else
            a=[1 5 9 13 17 21 25 29];
        end
        if ismember(ii,a)==0 %remove all of the interspersed S0, checks to see if index is in A, returns only false
            new_rf_wassr(jj)=rf_wassr(ii);
            CESTR_wassr(:,:,:,jj)=CEST_wassr(:,:,:,ii);
            jj=jj+1;
        else
            continue
        end
    end
    
    new_rf_wassr=new_rf_wassr.*128;
    if size(A_wassr,4) == 15
        xx_wassr=-0.5:.01:0.5;
    else
        xx_wassr=-1:.01:1;
    end
    
    xx_wassr_hz=xx_wassr.*128;
    A_wassr=double(A_wassr);
    
    % Spline Fit WASSR Data to find shift_w
    for slice=1:size(CESTR_wassr,3);
        for ii=1:size(CESTR_wassr,1);
            %disp(strcat('Processing Slice Number:  ',num2str(slice)))
            if mod(ii,20)==0
                disp(strcat('Fitting WASSR Row #: ', num2str(ii)))
            end
            for jj=1:size(CESTR_wassr,2);
                if maskCEST(ii,jj,slice)==1
                    dd=squeeze(double(CESTR_wassr(ii,jj,slice,:)));
                    ddSp=spline(new_rf_wassr, dd, xx_wassr_hz);
                    [c,I]=min(ddSp);
                    shift_w_wassr(ii, jj,slice,:)=xx_wassr_hz(I);
                else
                    shift_w_wassr(ii,jj,slice,:) = 0;
                end
            end
        end
    end
  
    
    c = clock;
    yy = num2str(c(1)); mm = num2str(c(2)); dd = num2str(c(3)); hh = num2str(c(4)); mmm = num2str(c(5));
    save_time = strcat(yy,mm,dd,'_',hh,mmm);
    
    
    savefile = strcat('WASSR_',save_time,'.mat');
    save(savefile, 'shift_w_wassr')

CESTR = single(CESTR);

%%


for slice=1:size(CESTR,3)
    
    
    disp(strcat('Shifting Slice Number:  ',num2str(slice)))
    
    for ii = 1:size(CESTR,1);
        
        if mod(ii, 20) == 0
            disp(strcat('Shifting Row Number:  ',num2str(ii)))
        end
        
        
        for jj = 1:size(CESTR,2);
            
            if maskCEST(ii,jj,slice) ~=0
                
                dd= double(squeeze(CESTR(ii,jj,slice,:)));
                
                rf_A_centered = rf_A - shift_w_wassr(ii,jj,slice);
                
                shifted_cest_data(ii,jj,slice,:) = interp1(rf_A_centered, dd, rf_A, 'pchip','extrap'); %extrap back to original x-values
                
                
                clear dd rr ff rf_A_centered
                
            else
                
                shifted_cest_data(ii,jj,slice,:) = zeros(length(rf_A),1);
                
                
            end
            
        end
    end
end

c = clock;
yy = num2str(c(1)); mm = num2str(c(2)); dd = num2str(c(3)); hh = num2str(c(4)); mmm = num2str(c(5));
save_time = strcat(yy,mm,dd,'_',hh,mmm);

savefile = strcat('Shifted_Fitted_CEST_Data_',save_time,'.mat');
save(savefile, 'shifted_cest_data', 'rf_A', 'shift_w_wassr', 'maskCEST')


%%
%%%Smooth shifted data%%%
disp('Smoothing z-spectra')
%keep shifted_cest_data shifted_residual_data rf_A mask shifted_fitted_cest_data CESTR
for slice = 1:size(CESTR,3);
    disp(strcat('Processing Slice Number:  ',num2str(slice)))
    
    for ii = 1:size(shifted_cest_data,1);
        
        if mod(ii, 20) == 0
            disp(strcat('Processing Row Number:  ',num2str(ii)))
        end
        
        for jj = 1:size(shifted_cest_data,2);
            
            if maskCEST(ii,jj, slice) == 1
                
                
                dd = double(squeeze(shifted_cest_data(ii,jj,slice,:)));
                if size(A,4) == 37
                    CESTR_smooth(ii,jj,slice,:) =[smooth(dd(1:5),10);dd(6:10);smooth(dd(11:end),5)]';
                else
                    CESTR_smooth(ii,jj,slice,:) =[smooth(dd(1:6),10);dd(7:15);smooth(dd(16:end),5)]';
                end
                %                CESTR_smooth(ii,jj,slice,:) =[smooth(dd(1:5),10);dd(6:10);smooth(dd(11:end),5)]'; %added 12/15/2015 sb
                %CESTR_smooth(ii,jj,slice,:) =[dd(1:6);dd(7:15);smooth(dd(16:end),5)]';
            else
                CESTR_smooth(ii,jj,slice,:) = zeros(1,size(shifted_cest_data,4));
                
            end
        end
    end
end


%%
Zspect = zeros(size(rf_A));
stdev = Zspect;

%Applymask Application
for ii = 1:size(CESTR_smooth, 4)
    [Zspect(ii), stdev(ii)] = applymask(squeeze(CESTR_smooth(:,:,1,ii)), maskCEST);
    
end
figure();
plot(rf_A./128, Zspect)
hold on
% legend('New', 'Old', 'Location', 'SE')
set(gca,'xdir','reverse')
title(sprintf('Z Spectrum for %s',ScanNumber))

c = clock;
yy = num2str(c(1)); mm = num2str(c(2)); dd = num2str(c(3)); hh = num2str(c(4)); mmm = num2str(c(5));
save_time = strcat(yy,mm,dd,'_',hh,mmm);

savefile = strcat('SmoothData_',save_time,'.mat');
save(savefile, 'CESTR_smooth','Zspect')

%% integrated MT asymmetry
if size(A,4) == 37
    integration_range=17:23; %(3.2-3.8)
else
    integration_range=26:32;
end

tt = zeros(size(CESTR_smooth,1),size(CESTR_smooth,2),size(CESTR_smooth,3));
for ii=1:size(CESTR_smooth,3)
    apt(:,:,ii) = trapz(rf_A(integration_range)./128, 1-CESTR_smooth(:,:,ii,integration_range),4); %1-CESTR_smooth
    mt(:,:,ii) = trapz(rf_A(3:5)./128, 1-CESTR_smooth(:,:,ii,3:5),4);
    tt(:,:,ii)=mt(:,:,ii)-apt(:,:,ii);
end

APTasym =tt;


MTRrex = zeros(size(CESTR_smooth,1),size(CESTR_smooth,2),size(CESTR_smooth,3));
for ii=1:size(CESTR_smooth,3)
    arex(:,:,ii) = 1./trapz(rf_A(integration_range)./128, CESTR_smooth(:,:,ii,integration_range),4); %1-CESTR_smooth
    mt1(:,:,ii) = 1./trapz(rf_A(3:5)./128, CESTR_smooth(:,:,ii,3:5),4);
    MTRrex(:,:,ii)=arex(:,:,ii)-mt1(:,:,ii);
end

MTRrex(isnan(MTRrex)) = 0;

c = clock;
yy = num2str(c(1)); mm = num2str(c(2)); dd = num2str(c(3)); hh = num2str(c(4)); mmm = num2str(c(5));
save_time = strcat(yy,mm,dd,'_',hh,mmm);

savefile = strcat('APT_',save_time,'.mat');
save(savefile, 'APTasym', 'integration_range', 'apt', 'mt','MTRrex');
