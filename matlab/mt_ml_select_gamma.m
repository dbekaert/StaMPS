function []=mt_ml_select_gamma(coh_thresh,image_fraction,weed_zero_elevation,endian,list)
%MT_ML_SELECT_GAMMA small baseline multilooked select for Gamma processed input
%  []=mt_ml_select_gamma(coh_thresh,image_fraction,weed_zero_elevation,list)
%
%   INPUTS: 
%   - coh_thresh: Threshold coherence. By default this value is 0.25
%   - image_fraction: fraction of images for a pixel that has met the
%                     coh_thresh value. Between 0-1, default this is 0.3
%   - weed_zero_elevation: default 'n' this is not done
%   - 
%   - list: txt file with ifgs folders that needs to be considered in the
%           processing. By default all date folders are considered. This is
%           a file given as DATE1_DATE2
%
%   This script is run in place of steps 1 to 5.
%
%   Andy Hooper, Oct 2015
%
% ===========================================================================
% Expected directory structure:
%
% SB:
%   SMALL_BASELINES/YYYYMMDD_YYYYMMDD/*.diff
%   SMALL_BASELINES/YYYYMMDD_YYYYMMDD/*.cc
%   SMALL_BASELINES/YYYYMMDD_YYYYMMDD/*.base
%   slc/*.mli  (single master)
%   slc/*.mli.par  (single master)
%   rslc/*.mli (single master slaves)
%   geo/*dem.rdc (1 file)
%   geo/*.lon (1 file)
%   geo/*.lat (1 file)
%   
%   Endianess should be the same for all bimnary files
%   Precision format should be the same for all binary files (e.g. FLOAT)
%   
% ===========================================================================

  
if nargin<1 || isempty(coh_thresh)
   coh_thresh=0.25;
end

if nargin<2 || isempty(image_fraction)
   image_fraction=0.3;
end

if nargin<3 || isempty(weed_zero_elevation)
    weed_zero_elevation='n';
end

if nargin<4
   endian='b';
end

if nargin<5 || isempty(list)   
    flag_all_ifgs = 1;
else
   flag_all_ifgs = 0;
end

%NB IFGS assumed in ascending date order ??
logit;

% Get matlab version as function arguments change with the matlab version
matlab_version = version('-release');           % [DB] getting the matlab version
matlab_version = str2num(matlab_version(1:4));  % [DB] the year


% Giving an output of the selected values:
logit(['coh_thresh : ' , num2str(coh_thresh)]);
logit(['image_fraction : ' , num2str(image_fraction)]);
logit(['weed_zero_elevation : ' , weed_zero_elevation]);


logit('Selecting multilooked pixels...')

PWD=pwd;
if ~strcmp(PWD(end-14:end),'SMALL_BASELINES')
    error('You are not in a SMALL_BASELINES directory');
else
    single_master_flag=0;
end

masterslc=dir(['..',filesep,'SLC',filesep,'*.mli']);
masterslc=masterslc(1).name;
master_date=masterslc(1:8);

fid=fopen('processor.txt','w');
fprintf(fid,'gamma');
fclose(fid);

rslcpar=dir(['..',filesep,'SLC',filesep,num2str(master_date),'*.par']);
rslcpar=['..',filesep,'SLC',filesep,rslcpar(1).name];

rps=readparm(rslcpar,'range_pixel_spacing');
imformat=readparm(rslcpar,'image_format');
rgn=readparm(rslcpar,'near_range_slc');
se=readparm(rslcpar,'sar_to_earth_center');
re=readparm(rslcpar,'earth_radius_below_sensor');
rgc=readparm(rslcpar,'center_range_slc');
n_laz=readparm(rslcpar,'azimuth_lines');
n_lrg=readparm(rslcpar,'range_samples');
prf=readparm(rslcpar,'prf');
looks=readparm(rslcpar,'range_looks');
looks_az=readparm(rslcpar,'azimuth_looks');
heading=readparm(rslcpar,'heading:');
sensor=readparm(rslcpar,'sensor');

cohname='.cc';
% allow a processing list of ifgs folders to be specified
if flag_all_ifgs==1
    dirs=ls(['*',filesep,'*',cohname]);
    dirs=strread(dirs,'%s');
    n_d=length(dirs); 
else
    dirs = textread(list,'%s');
    n_d = length(dirs); 
    for k=1:n_d
       dirs{k} = [dirs{k},filesep,cohname];
    end
end
    
ix=logical(zeros(n_laz,n_lrg,n_d));
i1=0;
rubbish=char([27,91,48,48,109]);

for i=1:n_d
    if strfind(dirs{i},cohname)
        i1=i1+1;
        filename=dirs{i};
        coh_ix=strfind(filename,rubbish);
        if ~isempty(coh_ix)
            filename=filename(coh_ix(end-1)+5:coh_ix(end)-1); % drop extra stuff sometimes added to start and end
        end
        fid=fopen(filename,'r',endian);
        if fid<=0
           error([filename,' does not appear to exist'])
        end
        logit(sprintf('reading %s',filename));
        if strcmpi(imformat,'FLOAT')
            coh=fread(fid,[n_lrg,inf],'float=>single');
        else
            error([imformat,' is not coded for'])
        end
        coh=coh';
        ix(:,:,i1)=coh>=coh_thresh;
        fclose(fid);
    end
end

ix=ix(:,:,1:i1);
ixsum=sum(ix,3);
ix=(ixsum>=i1*image_fraction);
[I,J]=find(ix);
ix=find(ix);
n_ps=length(ix);
logit(sprintf('%d pixels selected initially',n_ps))

ixsum=ixsum(ix);

[ny,nx]=size(coh);

lonname=dir(['..',filesep,'geo',filesep,'*.lon']);
if length(lonname)>0
    lonname=['..',filesep,'geo',filesep,lonname(1).name];
    logit(sprintf('reading %s',lonname));
    fid=fopen(lonname,'r',endian);
    %fid=fopen(lonname,'r');%%% FIX FOR NOW

    if strcmpi(imformat,'FLOAT')
        lon=fread(fid,[nx,ny],'float=>single');
    else
        error('unrecognised format')
    end
    lon=lon';
    fclose(fid);
else
    lon=repmat([1:n_lrg]*5/1e5,n_laz,1); % just for testing purposes
end
lon=lon(ix);

latname=dir(['..',filesep,'geo',filesep,'*.lat']);
if length(latname)>0
    latname=['..',filesep,'geo',filesep,latname(1).name];
    logit(sprintf('reading %s',latname));
    fid=fopen(latname,'r',endian);
    %fid=fopen(latname,'r'); %%% FIX FOR NOW

    if strcmpi(imformat,'FLOAT')
        lat=fread(fid,[nx,ny],'float=>single');
    else
        error('unrecognised format')
    end
    lat=lat';
    fclose(fid);
else
    lat=repmat([1:n_laz]'*20/1e5,1,n_lrg); % just for testing purposes
end
lat=lat(ix);

lonlat=[lon,lat];
clear lon lat

%%%Drop pixels with duplicate lon/lat coords
[dummy,ui]=unique(lonlat,'rows');
ix_weed=[1:n_ps]';
dups=setxor(ui,ix_weed); % pixels with duplicate lon/lat
ix_weed=logical(ones(n_ps,1));
for i=1:length(dups)
    dups_ix=find(lonlat(:,1)==lonlat(dups(i),1)&lonlat(:,2)==lonlat(dups(i),2));
    [dummy,max_i]=max(ixsum(dups_ix));
    ix_weed(dups_ix(setxor(max_i,[1:length(dups_ix)])))=0; % drop dups with lowest consistency
end
if ~isempty(dups)
    ix=ix(ix_weed);
    n_ps=length(ix);
    lonlat=lonlat(ix_weed,:);
    logit(sprintf('%d pixels with duplicate lon/lat dropped',length(dups)'));
end


ll0=(max(lonlat)+min(lonlat))/2;
xy=llh2local(lonlat',ll0)*1000;
xy=xy';
sort_x=sortrows(xy,1);
sort_y=sortrows(xy,2);
n_pc=round(n_ps*0.001);
bl=mean(sort_x(1:n_pc,:)); % bottom left corner
tr=mean(sort_x(end-n_pc:end,:)); % top right corner
br=mean(sort_y(1:n_pc,:)); % bottom right  corner
tl=mean(sort_y(end-n_pc:end,:)); % top left corner

theta=(180-heading)*pi/180;
if theta>pi
    theta=theta-2*pi;
end

rotm=[cos(theta),sin(theta); -sin(theta),cos(theta)];
xy=xy';
xynew=rotm*xy; % rotate so that scene axes approx align with x=0 and y=0
if max(xynew(1,:))-min(xynew(1,:))<max(xy(1,:))-min(xy(1,:)) &...
   max(xynew(2,:))-min(xynew(2,:))<max(xy(2,:))-min(xy(2,:))
    xy=xynew; % check that rotation is an improvement
    logit(['Rotating by ',num2str(theta*180/pi),' degrees']);
end
        
xy=single(xy');
[dummy,sort_ix]=sortrows(xy,[2,1]); % sort in ascending y order
xy=xy(sort_ix,:);
xy=[[1:n_ps]',xy];
xy(:,2:3)=round(xy(:,2:3)*1000)/1000; % round to mm

lonlat=lonlat(sort_ix,:);
ix=ix(sort_ix);
ij=[[1:n_ps]',I(sort_ix)*looks_az-round(looks_az/2),J(sort_ix)*looks-round(looks/2)];



if flag_all_ifgs==1
    dirs=ls('-1',['*',filesep,'*.diff']);
    dirs=strread(dirs,'%s');
    n_ifg=length(dirs);
else
    dirs = textread(list,'%s');
    n_ifg = length(dirs); 
    for k=1:n_ifg
       dirs{k} = [dirs{k} filesep dirs{k},'.diff'];
    end
end

ph=zeros(n_ps,n_ifg,'single');
coh_ps=zeros(n_ps,n_ifg,'single');
day1_yyyymmdd=zeros(n_ifg,1);
day2_yyyymmdd=zeros(n_ifg,1);
rg=rgn+ij(:,3)*rps;
look=acos((se^2+rg.^2-re^2)./(2*se*rg)); % Satellite look angles
inci=acos((se^2-re^2-rg.^2)./(2*re*rg));
mean_incidence=mean(inci);
mean_range=rgc;
mean_az=n_laz/2-0.5; 
bperp_mat=zeros(n_ps,n_ifg,'single');

i1=0;
for i=1:n_ifg
    if strfind(dirs{i},'.diff')
        i1=i1+1;
        filename=dirs{i};
        cint_ix=strfind(filename,rubbish);
        if ~isempty(cint_ix)
            filename=filename(cint_ix(end-1)+5:cint_ix(end)-1);
        end
        if single_master_flag==0
            day1_yyyymmdd(i1)=str2num(filename(1:8));
            day2_yyyymmdd(i1)=str2num(filename(10:17));
        else
            day1_yyyymmdd(i1)=str2num(masterdate);
            day2_yyyymmdd(i1)=str2num(filename(1:8));
        end
        
        logit(sprintf('reading %s',filename));
        if strcmpi(imformat,'FLOAT')
            ifg=readcpx(filename,n_lrg,n_laz,endian);
        else
            error('unrecognised format')
        end
            
        ph(:,i1)=ifg(ix);
        
        filename=[filename(1:end-5),cohname];
        logit(sprintf('reading %s',filename));
        fid=fopen(filename,'r',endian);
        if strcmpi(imformat,'FLOAT')
            coh=fread(fid,[nx,ny],'float=>single');
        else
            error('unrecognised format')
        end
        coh=coh';
        fclose(fid);
        coh_ps(:,i1)=coh(ix);
        
        basename=[filename(1:end-3),'.base'];
        B_TCN=readparm(basename,'initial_baseline(TCN):',3,0);
        BR_TCN=readparm(basename,'initial_baseline_rate:',3,0);
        bc=B_TCN{2}+BR_TCN{2}*(ij(:,2)-mean_az)/prf;
        bn=B_TCN{3}+BR_TCN{3}*(ij(:,2)-mean_az)/prf;
        bperp_mat(:,i)=bc.*cos(look)-bn.*sin(look); % Convert baselines from (T)CN to perp-para coordinates
    end
end

bperp=mean(bperp_mat)';
if single_master_flag~=0
    bperp_mat=bperp_mat(:,[1:master_ix-1,master_ix+1:end]); % depends on whether null interferogram included in list?
end

n_ifg=i1;
ph=ph(:,1:i1);
day1_yyyymmdd=day1_yyyymmdd(1:i1);
day2_yyyymmdd=day2_yyyymmdd(1:i1);

nzix=sum(ph==0,2)==0;
ph=ph(nzix,:);
xy=xy(nzix,:);
ij=ij(nzix,:);
lonlat=lonlat(nzix,:);
bperp_mat=bperp_mat(nzix,:);
ix=ix(nzix);
n_ps=size(ph,1);


year=floor(day1_yyyymmdd/10000);
month=floor((day1_yyyymmdd-year*10000)/100);
monthday=day1_yyyymmdd-year*10000-month*100;
ifgday(:,1)=datenum(year,month,monthday);
year=floor(day2_yyyymmdd/10000);
month=floor((day2_yyyymmdd-year*10000)/100);
monthday=day2_yyyymmdd-year*10000-month*100;
ifgday(:,2)=datenum(year,month,monthday);
[day,I,J]=unique(ifgday(:));
ifgday_ix=reshape(J,n_ifg,2);
master_day_yyyymmdd=str2num(master_date);
year=floor(master_day_yyyymmdd/10000);
month=floor((master_day_yyyymmdd-year*10000)/100);
monthday=master_day_yyyymmdd-year*10000-month*100;
master_day=datenum(year,month,monthday);
master_ix=sum(master_day>day)+1;
n_image=size(day,1);
psver=2;
save psver psver

demname=dir(['..',filesep,'geo',filesep,'*dem.rdc']);
demname=['..',filesep,'geo',filesep,demname(1).name];
hgt=zeros(ny,nx,'single');


logit(sprintf('reading %s',demname));
fid=fopen(demname,'r',endian);
if strcmpi(imformat,'FLOAT')
    hgt=fread(fid,[nx,ny],'float=>single');
else
    error('unrecognised format')
end
hgt=hgt';

fclose(fid);
hgt=hgt(ix);


% weeding zero height information when requested
if ~strcmpi(weed_zero_elevation,'n')
    nzix=hgt>0;
    hgt=hgt(nzix);
    ph=ph(nzix,:);
    coh_ps= coh_ps(nzix,:);
    bperp_mat=bperp_mat(nzix,:);
    xy=xy(nzix,:);
    ij=ij(nzix,:);
    lonlat=lonlat(nzix,:);
    ix=ix(nzix);
    n_ps=size(ph,1);
end
save('hgt2','hgt')
save('bp2','bperp_mat')




if single_master_flag~=0
    ph_rc=[ph(:,1:master_ix-1),zeros(size(ph,1),1),ph(:,master_ix:end)];
    n_ifg=n_ifg+1;
    coh_ps= [coh_ps(:,1:master_ix-1),zeros(size(coh_ps,1),1),coh_ps(:,master_ix:end)];
else
    ph_rc=ph;
end
save('rc2','ph_rc');
save('pm2','coh_ps')

if single_master_flag==0
    sb_parms_initial
else
    ps_parms_initial
end


setparm('heading',heading,1);
setparm('weed_zero_elevation',weed_zero_elevation)

freq=readparm(rslcpar,'radar_frequency:');
lambda=299792458/freq;
setparm('lambda',lambda,1);

if ~isempty(strfind('sensor','ASAR'))
    platform='ENVISAT';
else
    platform=sensor; % S1A for Sentinel-1A
end
setparm('platform',platform,1);



% constructing ps2.mat
save('ps2','ij','lonlat','xy','day','ifgday','ifgday_ix','bperp','master_day','master_ix','n_ifg','n_image','n_ps','ll0','master_ix','mean_incidence','mean_range');


logit(1);
