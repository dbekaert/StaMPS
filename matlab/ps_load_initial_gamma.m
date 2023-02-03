function []=ps_load_initial_gamma(endian)
%PS_LOAD_INITIAL_GAMMA Initial load of files into matlab workspace
%
%   Andy Hooper, Dec 2012
%
%   =======================================================================
%   01/2016 DB: Replace save with stamps_save which checks for var size when
%               saving 
%   01/2016 AH: Save sensor
%   12/2016 AH: Set xy from lonlat 
%   08/2017 AH: Allow for missing null ifg
%   =======================================================================
  

if nargin<1
   endian='b';
end

%NB IFGS assumed in ascending date order
logit;
logit('Loading data into matlab...')

phname=['./pscands.1.ph']; % for each PS candidate, a float complex value for each ifg
ijname=['./pscands.1.ij']; % ID# Azimuth# Range# 1 line per PS candidate
llname=['./pscands.1.ll']; % output from pt2geo
xyname=['./pscands.1.xy']; % output from pt2geo
hgtname=['./pscands.1.hgt']; % 
daname=['./pscands.1.da']; % 
rscname = ['../rsc.txt'];
pscname=['../pscphase.in'];


psver=1;
fid=fopen(rscname);
if fid<0
    error([rscname,' does not exist'])
end
rslcpar=textscan(fid,'%s');
rslcpar=rslcpar{1}{1};
fclose(fid);

fid=fopen(pscname);
if fid<0
    error([pscname,' does not exist'])
end
ifgs=textscan(fid,'%s');
fclose(fid);
ifgs=ifgs{1}(2:end);
nb=length(ifgs{1});
master_day=str2num(ifgs{1}(nb-27:nb-20));
n_ifg=length(ifgs);
n_image=n_ifg;
day=zeros(n_ifg,1);
for i=1:n_ifg
    day(i)=str2num(ifgs{i}(nb-12:nb-5));
end

year=floor(day/10000);
month=floor((day-year*10000)/100);
monthday=day-year*10000-month*100;
day=datenum(year,month,monthday);

master_day_yyyymmdd=master_day;
year=floor(master_day/10000);
month=floor((master_day-year*10000)/100);
monthday=master_day-year*10000-month*100;
master_day=datenum(year,month,monthday);

master_ix=sum(day<master_day)+1;
if day(master_ix)~=master_day
    master_master_flag='0' % no null master-master ifg provided
    day=[day(1:master_ix-1);master_day;day(master_ix:end)];
else
    master_master_flag='1' % yes, null master-master ifg provided
end

heading=readparm(rslcpar,'heading:');
setparm('heading',heading,1);

freq=readparm(rslcpar,'radar_frequency:');
lambda=299792458/freq;
setparm('lambda',lambda,1);

sensor=readparm(rslcpar,'sensor:');
if ~isempty(strfind('sensor','ASAR'))
    platform='ENVISAT';
else
    platform=sensor; % S1A for Sentinel-1A
end
setparm('platform',platform,1);



ij=load(ijname);
n_ps=size(ij,1);

rps=readparm(rslcpar,'range_pixel_spacing');
rgn=readparm(rslcpar,'near_range_slc');
se=readparm(rslcpar,'sar_to_earth_center');
re=readparm(rslcpar,'earth_radius_below_sensor');
rgc=readparm(rslcpar,'center_range_slc');
naz=readparm(rslcpar,'azimuth_lines');
prf=readparm(rslcpar,'prf');


mean_az=naz/2-0.5; % mean azimuth line

rg=rgn+ij(:,3)*rps;
look=acos((se^2+rg.^2-re^2)./(2*se*rg)); % Satellite look angles

bperp_mat=zeros(n_ps,n_image,'single');
for i=1:n_ifg
    basename=[ifgs{i}(1:nb-4),'base'];
    B_TCN=readparm(basename,'initial_baseline(TCN):',3);
    BR_TCN=readparm(basename,'initial_baseline_rate:',3);
    bc=B_TCN{2}+BR_TCN{2}*(ij(:,2)-mean_az)/prf;
    bn=B_TCN{3}+BR_TCN{3}*(ij(:,2)-mean_az)/prf;
    bperp_mat(:,i)=bc.*cos(look)-bn.*sin(look); % Convert baselines from (T)CN to perp-para coordinates
    %bpara=bc*sin(look)+bn*cos(look)
end

bperp=mean(bperp_mat)';
if master_master_flag==1
    bperp_mat=bperp_mat(:,[1:master_ix-1,master_ix+1:end]);
else
    bperp=[bperp(1:master_ix-1);0;bperp(master_ix:end)];
end
%bperp=[bperp(1:master_ix-1);0;bperp(master_ix:end)]; % insert master-master bperp (zero)
%bperp_mat=repmat(single(bperp([1:master_ix-1,master_ix+1:end]))',n_ps,1);

inci=acos((se^2-re^2-rg.^2)./(2*re*rg));
mean_incidence=mean(inci);
mean_range=rgc;
    
fid=fopen(phname,'r',endian);
ph=zeros(n_ps,n_ifg,'single');
byte_count=n_ps*2;
for i=1:n_ifg
    [ph_bit,byte_count]=fread(fid,[(n_ps)*2,1],'float=>single');
    ph(:,i)=complex(ph_bit(1:2:end),ph_bit(2:2:end));
end

zero_ph=sum(ph==0,2);
nonzero_ix=zero_ph<=1;       % if more than 1 phase is zero, drop node
if master_master_flag==1
    ph(:,master_ix)=1;
else
    ph=[ph(:,1:master_ix-1),ones(n_ps,1),ph(:,master_ix:end)];
    n_ifg=n_ifg+1;
    n_image=n_image+1;
end

if exist(llname,'file')
    fid=fopen(llname,'r');
    lonlat=fread(fid,[2,inf],'float',endian);
    lonlat=lonlat';
    fclose(fid);
else
    error([llname,' does not exist']);
end

ll0=(max(lonlat)+min(lonlat))/2;

xy=llh2local(lonlat',ll0)'*1000;

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
    disp(['Rotating by ',num2str(theta*180/pi),' degrees']);
end
        
xy=single(xy');
[dummy,sort_ix]=sortrows(xy,[2,1]); % sort in ascending y order
xy=xy(sort_ix,:);
xy=[[1:n_ps]',xy];
xy(:,2:3)=round(xy(:,2:3)*1000)/1000; % round to mm

ph=ph(sort_ix,:);
ij=ij(sort_ix,:);
ij(:,1)=1:n_ps;
lonlat=lonlat(sort_ix,:);
bperp_mat=bperp_mat(sort_ix,:);


% removing Nan values if pressent in lonlat, phase
ix_nan = (sum(isnan(lonlat),2)>=1 | sum(isnan(ph),2)>=1);
lonlat(ix_nan,:)=[];
ij(ix_nan,:)=[];
xy(ix_nan,:)=[];
n_ps = size(lonlat,1);
% update the first column
ij(:,1)=[1:n_ps]';
xy(:,1)=[1:n_ps]';

savename=['ps',num2str(psver)];
stamps_save(savename,ij,lonlat,xy,bperp,day,master_day,master_ix,n_ifg,n_image,n_ps,sort_ix,ll0,master_ix,mean_incidence,mean_range);

save psver psver

phsavename=['ph',num2str(psver)];
ph(ix_nan,:)=[];
stamps_save(phsavename,ph);

bpsavename=['bp',num2str(psver)];
bperp_mat(ix_nan,:)=[];
stamps_save(bpsavename,bperp_mat);

lasavename=['la',num2str(psver)];
la=inci(sort_ix); % store incidence not look angle for gamma
la(ix_nan,:)=[];
stamps_save(lasavename,la);

if exist(daname,'file')
  D_A=load(daname);
  D_A=D_A(sort_ix);
  D_A(ix_nan,:)=[];
  dasavename=['da',num2str(psver)];
  stamps_save(dasavename,D_A);
end

if exist(hgtname,'file')
    fid=fopen(hgtname,'r');
    hgt=fread(fid,[1,inf],'float',endian);
    hgt=hgt';
    hgt=hgt(sort_ix);
    fclose(fid);
    hgt(ix_nan,:)=[];
    hgtsavename=['hgt',num2str(psver)];
    stamps_save(hgtsavename,hgt);
end

logit(1);
