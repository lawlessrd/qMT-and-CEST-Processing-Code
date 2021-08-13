function [rvt, actual_respWASSR]=getRespRVT_r55_3TA(wassr_check);
%04/01/2016
%%  Read the respiratory data
%   This is for r5.3, where the physlog data is slightly modified. 
%   Modified from getRespRVT, and originally from Rob B.'s vuRETROICOR function 
%   Reads in the scanner's log file. 
%   Normalizes the respiration bellows data on a scale of 0 to 1
%   KPO modified 6/26/19 because loadPARREC cannot open PARRECs from R5.5/3T-A
if strcmp(wassr_check, 'y')==1
    logpath=sprintf('%s/', pwd);
    c=dir('SCAN*_WASSR.log');
    logname=c.name;%'SCANPHYSLOG20160109134009_CEST.log';
    d=dir('Smith*WASSR*.PAR');
    parpath=sprintf('%s/', pwd);
    parname=d.name; %'Smith_227055_WIP_150ms_5NSA_SENSE_11_1.PAR'; %need to get the parameter inputs
    parfile=[parpath parname];
    logfile=[logpath logname];

    info = loadParRec(parfile);
    num_slices = info.datasize(1,3);
    scant = info.imgdef.dyn_scan_begin_time.vals+0.075;
    if size(scant,2) > size(scant,1), scant = scant'; end
    vat = scant(end,1) - scant(end-1,1);
    if vat == 0
        vat = scant(end,1) - scant(end-num_slices,1);
    end
    % added 2015-07-28 to fix incorrect vat calculation
    if vat == 0
        vat = scant(end,1) - scant(end-num_slices*size(info.imgdef.image_type_mr.uniq,1),1);
    end

    dynamics = info.Parms.max_dynamics;

    resp = phys_filt_sync_ricor_r53(logfile, 6, vat, dynamics);
    %card = phys_filt_sync_ricor(logfile, 5, vat, dynamics);
    resp_norm=(resp-min(resp))/range(resp);

    % 01/07/2016
    acquisition_time=info.Parms.scan_dur; %seconds %166.7;%84.7 for wassr, 824 for cest
    resp_sampling=size(resp,2);
    x=1:resp_sampling;
    sr=acquisition_time/resp_sampling;
    time=x.*sr;

    SamplesPerSec=size(time,2)/acquisition_time;


     for ii=1:dynamics
         a=scant(ii)*SamplesPerSec;
         actual_respWASSR(ii)=resp_norm(ceil(a));
     end


keep actual_respWASSR wassr_check
else
    disp('No WASSR Resp');
end
 %%
disp('Getting CEST data');
logpath=sprintf('%s/', pwd);
c=dir('SCAN*_CEST.log');
logname=c.name;
%logname='SCANPHYSLOG20160330124140_CESTsc12.log';%'SCANPHYSLOG20160109134009_CEST.log'; %c.name;

d=dir('Smith*150ms_*NSA*.PAR'); %KPO removed 5 from '5NSA' to processed 3D CEST w/1NSA (same TR) 6/26/19
parpath=sprintf('%s/', pwd);
%parname='Smith_20160330_egg_WIP_150ms_5NSA_2uT_SENSE_12_1.PAR';
parname=d.name;%'Smith_227055_WIP_150ms_5NSA_SENSE_11_1.PAR'; %need to get the parameter inputs
parfile=[parpath parname];
logfile=[logpath logname];

% info = loadParRec(parfile); %KPO commented out 6/26/19 to rewrite for R5.5
% num_slices = info.datasize(1,3);
% scant = info.imgdef.dyn_scan_begin_time.vals +.305; % 0.305 is our TR in seconds

%KPO modified to use parrec from 3T-A in R5.5
info = vuOpenImage(parfile);
num_slices = info.Parms.max_slices;
slice_start_t = info.Parms.tags(:,32);
dyns = 1:num_slices:(size(slice_start_t,1)-(num_slices-1));
TR = info.Parms.repetition_time/1000; %TR in seconds
scant = slice_start_t(dyns) + TR;


if size(scant,2) > size(scant,1), scant = scant'; end
vat = scant(end,1) - scant(end-1,1);
if vat == 0
    vat = scant(end,1) - scant(end-num_slices,1);
end
% added 2015-07-28 to fix incorrect vat calculation
if vat == 0
    vat = scant(end,1) - scant(end-num_slices*size(info.imgdef.image_type_mr.uniq,1),1); %KPO need to edit for R5.5
end

dynamics = info.Parms.max_dynamics;

resp = phys_filt_sync_ricor_r53(logfile, 6, vat, dynamics);
%card = phys_filt_sync_ricor(logfile, 5, vat, dynamics);
resp_norm=(resp-min(resp))/range(resp);


acquisition_time=info.Parms.scan_dur;%166.7;%84.7 for wassr, 824 for cest
resp_sampling=size(resp,2);
x=1:resp_sampling;
sr=acquisition_time/resp_sampling;
time=x.*sr;

SamplesPerSec=size(time,2)/acquisition_time;


NumPE= info.Parms.scan_resolution(2);
EPIfactor= info.Parms.epi_factor;
SENSEFactor= 2;
NSA=info.Parms.tags(1,35);
    
linesPerDynamic=(NumPE/EPIfactor/SENSEFactor)*NSA; %NumPE/EPIFactor/SENSEFactor * NSA
numShots=(NumPE/EPIfactor/SENSEFactor);
kk=1;
for jj=1:dynamics
    for ii=1:linesPerDynamic %5 averages %%changed here
        new(kk)=scant(jj)+TR*(ii-1);
        kk=kk+1;
    end
end

for ii=1:size(new,2)
     a=new(ii)*SamplesPerSec; % %for this, scant has to be +TR
     actual_resp(ii)=resp_norm(ceil(a));
end



% % averages of shots are acquired right after one another, RVT: see Birn et al. 2006
for ii=1:dynamics

   kk=1;

    [pks, locs]=findpeaks(actual_resp(((ii-1)*linesPerDynamic) +1:ii*linesPerDynamic)); 
    numCycles=size(pks,2);
    jj=1;
    for peaks=1:size(locs,2)-1;
        dt(jj)=new(locs(peaks+1))-new(locs(peaks));
        jj=jj+1;
    end
    range_data(ii)=max(actual_resp(((ii-1)*linesPerDynamic )+1:ii*linesPerDynamic))-min(actual_resp(((ii-1)*linesPerDynamic )+1:ii*linesPerDynamic)); %./mean(actual_resp(((ii-1)*linesPerDynamic )+1:ii*linesPerDynamic));

    rvt(ii)=(range_data(ii)./mean(dt)); %range over all the shots in one dynamic/ mean time between peaks

end

if strcmp(wassr_check, 'y')==1
    keep rvt actual_respWASSR
else 
    keep rvt
end
