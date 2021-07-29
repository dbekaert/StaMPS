function []=ps_load_initial_gsar(endian)
%PS_LOAD_INITIAL_GSAR Initial load of files into matlab workspace
%
%   Andy Hooper, Sep 2015
%
%   =======================================================================
%   =======================================================================
  

if nargin<1
   endian='l';
end

%NB IFGS assumed in ascending date order
logit;
logit('Loading data into matlab...')

phname=['./pscands.1.ph']; % for each PS candidate, a float complex value for each ifg
ijname=['./pscands.1.ij']; % ID# Azimuth# Range# 1 line per PS candidate
lonname=['./pscands.1.lon']; % 
latname=['./pscands.1.lat']; % 
hgtname=['./pscands.1.hgt']; % 
incname=['./pscands.1.inc']; % 
daname=['./pscands.1.da']; % 
rscname = ['../rsc.txt'];
pscname=['../pscphase.in'];
calname=['calamp.out'];             % amplitide calibrations

psver=1;

fid=fopen(pscname);
ifgs=textscan(fid,'%s');
fclose(fid);
ifgs=ifgs{1}(2:end);
nb=length(ifgs{1});
n_ifg=length(ifgs)+1;
n_image=n_ifg;
day1=zeros(n_ifg-1,1);
day2=day1;
bperp=zeros(n_ifg-1,1);
for i=1:n_ifg-1
    day2(i)=str2num(ifgs{1}(nb-19:nb-12));
    day1(i)=str2num(ifgs{i}(nb-35:nb-28));
    ifgxml=xml2struct([ifgs{i},'.xml']);
    [~,xmlchild,~]=ifgxml.Children.Children;
    names={xmlchild.Name};
    ix=find(strcmpi(names,'baseline'));
    bperp(i)=str2num(xmlchild(ix).Children.Data);
end
if sum(abs(diff(day1)))==0
    master_day=day1(1);
    day=day2;
    switch_sign_flag=0;
elseif sum(abs(diff(day2)))==0
    master_day=day2(1);
    day=day1;
    switch_sign_flag=1;
else
    logit('error: Interferograms are not single master')
    error('abending...')
end
    
year=floor(day/10000);
month=floor((day-year*10000)/100);
monthday=day-year*10000-month*100;
day=datenum(year,month,monthday);
[day,day_ix]=sort(day);
bperp=bperp(day_ix);

master_day_yyyymmdd=master_day;
year=floor(master_day/10000);
month=floor((master_day-year*10000)/100);
monthday=master_day-year*10000-month*100;
master_day=datenum(year,month,monthday);

master_ix=sum(day<master_day)+1;
day=[day(1:master_ix-1);master_day;day(master_ix:end)]; % insert master day

ix=find(strcmpi(names,'trackangle'));
heading=str2num(xmlchild(ix).Children.Data); % ground track not satellite heading stored for  gsar
setparm('heading',heading,1);

ix=find(strcmpi(names,'radarfreq'));
freq=str2num(xmlchild(ix).Children.Data); % numeric
lambda=299792458/freq;
setparm('lambda',lambda,1);

ix=find(strcmpi(names,'platform'));
platform=xmlchild(ix).Children.Data; % string
setparm('platform',platform,1);

ij=load(ijname);
n_ps=size(ij,1);

bperp_mat=repmat(single(bperp)',n_ps,1); % create bperp matrix without master
bperp=[bperp(1:master_ix-1);0;bperp(master_ix:end)]; % insert master in bperp vector

if exist(calname,'file')
    [calfile,calconst]=textread(calname,'%s%f');
    caldate=zeros(length(calfile),1);
    for i = 1 : length(calfile)
        aa=strread(calfile{i},'%s','delimiter','/');
        bb=str2num(aa{end}(1:8));
        if isempty(bb)
            bb=str2num(aa{end-1}(1:8));
        end
        caldate(i)=bb;    end
    not_master_ix=caldate~=master_day_yyyymmdd;
    caldate=caldate(not_master_ix);
    calconst=calconst(not_master_ix);
    [Y,I]=sort(caldate);
    calconst=calconst(I);
else
    calconst=ones(n_ifg-1,1);
end


fid=fopen(incname,'r');
inc=fread(fid,[1,inf],'double',endian)'*pi/180; 
fclose(fid);
mean_incidence=mean(inc);
%mean_range=rgc;
mean_range=[];
    
fid=fopen(phname,'r',endian);
ph=zeros(n_ps,n_ifg,'single');
byte_count=n_ps*2;
for i=1:n_ifg-1
    [ph_bit,byte_count]=fread(fid,[(n_ps)*2,1],'float=>single');
    ph(:,i)=complex(ph_bit(1:2:end),ph_bit(2:2:end));
end

ph=ph(:,day_ix);
zero_ph=sum(ph==0,2);
nonzero_ix=zero_ph<=1;       % if more than 1 phase is zero, drop node
ph=ph./repmat(calconst',n_ps,1);  % scale amplitudes
if switch_sign_flag
    ph=conj(ph);
end
ph=[ph(:,1:master_ix-1),ones(n_ps,1),ph(:,master_ix:end)]; % insert zero phase master-master ifg

fid=fopen(lonname,'r');
lon=fread(fid,[1,inf],'double',endian);
lon=lon';
fclose(fid);

fid=fopen(latname,'r');
lat=fread(fid,[1,inf],'double',endian);
lat=lat';
fclose(fid);

lonlat=[lon,lat];

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


savename=['ps',num2str(psver)];

save(savename,'ij','lonlat','xy','bperp','day','master_day','master_ix','n_ifg','n_image','n_ps','sort_ix','ll0','master_ix','mean_incidence','mean_range');
save psver psver

phsavename=['ph',num2str(psver)];
save(phsavename,'ph');

bpsavename=['bp',num2str(psver)];
save(bpsavename,'bperp_mat');

lasavename=['la',num2str(psver)];
la=inc(sort_ix); % store incidence not look angle for GSAR
save(lasavename,'la');

if exist(daname,'file')
  D_A=load(daname);
  D_A=D_A(sort_ix);
  dasavename=['da',num2str(psver)];
  save(dasavename,'D_A');
end

if exist(hgtname,'file')
    fid=fopen(hgtname,'r');
    hgt=fread(fid,[1,inf],'float',endian);
    hgt=hgt';
    hgt=hgt(sort_ix);
    fclose(fid);
    hgtsavename=['hgt',num2str(psver)];
    save(hgtsavename,'hgt');
end

logit(1);
