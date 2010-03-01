function []=sb_load_initial()
%SB_LOAD_INITIAL Initial load of PS files into matlab workspace
%   SB_LOAD_INITIAL() loads 'pscands.1.*' files and stores them
%   in various workspaces. The version number, PSVER, is set to 1.
%
%   Andy Hooper, Spetember 2006
%
%   Change Log: 
%   11/2007 AH: change bperp to be bperp at 0 m

%NB IFGS assumed in ascending date order

phname=['pscands.1.ph'];            % for each PS candidate, a float complex value for each ifg
ijname=['pscands.1.ij'];            % ID# Azimuth# Range# 1 line per PS candidate
bperpname=['bperp.1.in'];           % in meters 1 line per slave image
dayname=['day.1.in'];               % YYYYMMDD, 1 line per slave image
ifgdayname=['ifgday.1.in'];         % YYYYMMDD YYYYMMDD, 1 line per ifg
masterdayname=['master_day.1.in'];  % YYYYMMDD
llname=['pscands.1.ll'];            % 2 float values (lon and lat) per PS candidate
daname=['pscands.1.da'];            % 1 float value per PS candidate
hgtname=['pscands.1.hgt'];          % 1 float value per PS candidate
laname=['look_angle.1.in'];         % grid of look angle values
headingname=['heading.1.in'];       % satellite heading
lambdaname=['lambda.1.in'];         % wavelength
calname=['calamp.out'];             % amplitide calibrations
widthname=['width.txt'];            % width of interferograms
lenname=['len.txt'];                % length of interferograms

psver=1;

if ~exist(dayname,'file')
    dayname= ['../',dayname];
end
day_yyyymmdd=load(dayname);
year=floor(day_yyyymmdd/10000);
month=floor((day_yyyymmdd-year*10000)/100);
monthday=day_yyyymmdd-year*10000-month*100;
slave_day=datenum(year,month,monthday);
[slave_day,day_ix]=sort(slave_day);
day_yyyymmdd=day_yyyymmdd(day_ix);


if ~exist(masterdayname,'file')
    masterdayname= ['../',masterdayname];
end
master_day_yyyymmdd=load(masterdayname);
year=floor(master_day_yyyymmdd/10000);
month=floor((master_day_yyyymmdd-year*10000)/100);
monthday=master_day_yyyymmdd-year*10000-month*100;
master_day=datenum(year,month,monthday);
master_ix=sum(master_day>slave_day)+1;
day=[slave_day(1:master_ix-1);master_day;slave_day(master_ix:end)]; % insert master day 
n_image=size(day,1);

if ~exist(ifgdayname,'file')
    ifgdayname= ['../',ifgdayname];
end
ifgday=load(ifgdayname);
year=floor(ifgday/10000);
month=floor((ifgday-year*10000)/100);
monthday=ifgday-year*10000-month*100;
ifgday=datenum(year,month,monthday);
n_ifg=size(ifgday,1);
[found_true,ifgday_ix]=ismember(ifgday,day);
if sum(found_true~=1)>0
   error('one or more days in ifgday.1.in not in day.1.in')
end

if ~exist(bperpname,'file')
    bperpname= ['../',bperpname];
end
bperp=load(bperpname); 
bperp=bperp(day_ix);
bperp=[bperp(1:master_ix-1);0;bperp(master_ix:end)]; % insert master-master bperp (zero)
bperp=bperp(ifgday_ix(:,2))-bperp(ifgday_ix(:,1));

 
if ~exist(headingname,'file')
    headingname= ['../',headingname];
end
heading=load(headingname);
if isempty(heading)
    error('heading.1.in is empty')
end
setparm('heading',heading,1);

if ~exist(lambdaname,'file')
    lambdaname= ['../',lambdaname];
end
lambda=load(lambdaname);
setparm('lambda',lambda,1);


ij=load(ijname);
n_ps=size(ij,1);

fid=fopen(phname,'r');
if fid < 0
   error([phname,' cannot be openned'])
end
ph_bit=fread(fid,[1,1],'float');
float_bytes=ftell(fid);
fseek(fid,0,1);
nbytes=ftell(fid);
if nbytes~=n_ps*n_ifg*float_bytes*2
    error([phname,' has ',num2str(nbytes/float_bytes/2),' complex float values which does not equal ',num2str(n_ps),'*',num2str(n_ifg),' (pixels*images)'])
end
fseek(fid,0,-1);
ph=zeros(n_ps,n_ifg,'single');
byte_count=n_ps*2;
for i=1:n_ifg
    [ph_bit,byte_count]=fread(fid,[(n_ps)*2,1],'float');
    ph_bit=single(ph_bit);
    ph(:,i)=complex(ph_bit(1:2:end),ph_bit(2:2:end));
end

fclose(fid);
%ph=single(ph');
%ph=complex(ph(:,1:2:end),ph(:,2:2:end));
zero_ph=sum(ph==0,2);
nonzero_ix=zero_ph<=1;       % if more than 1 phase is zero, drop node
ph(ph~=0)=ph(ph~=0)./abs(ph(ph~=0));  % scale amplitudes


fid=fopen(llname,'r');
if fid < 0
   error([llname,' cannot be openned'])
end
lonlat=fread(fid,[2,inf],'float');
lonlat=lonlat';
fclose(fid);
if size(lonlat,1) ~= n_ps
   error([llname,' has ',num2str(size(lonlat,1)),' entries whereas ',ijname,' has ',num2str(n_ps)])
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

savename=['ps',num2str(psver)];
phname=['ph',num2str(psver)];

save(savename,'ij','lonlat','xy','bperp','day','master_day','master_ix','ifgday','ifgday_ix','n_image','n_ifg','n_ps','sort_ix','ll0','master_ix','day_ix');
save(phname,'ph');
save psver psver
    
if exist(daname,'file')
  D_A=load(daname);
  if size(D_A,1) ~= n_ps
     error([daname,' has ',num2str(size(D_A,1)),' entries whereas ',ijname,' has ',num2str(n_ps)])
  end
  D_A=D_A(sort_ix);
  dasavename=['da',num2str(psver)];
  save(dasavename,'D_A');
end


if exist(hgtname,'file')
    fid=fopen(hgtname,'r');
    if fid < 0
       error([hgtname,' cannot be openned'])
    end
    hgt=fread(fid,[1,inf],'float');
    hgt=hgt';
    fclose(fid);
    if size(hgt,1) ~= n_ps
       error([hgtname,' has ',num2str(size(hgt,1)),' entries whereas ',ijname,' has ',num2str(n_ps)])
    end
    hgt=hgt(sort_ix);
    hgtsavename=['hgt',num2str(psver)];
    save(hgtsavename,'hgt');
end

if ~exist(widthname,'file')
    widthname= ['../',widthname];
end
width=load(widthname);

if ~exist(lenname,'file')
    lenname= ['../',lenname];
end
len=load(lenname);

[gridX,gridY]=meshgrid(linspace(0,width,50),linspace(0,len,50));

if ~exist(laname,'file')
    laname= ['../',laname];
end
if exist(laname,'file')
    look_angle=load(laname);
    la0=look_angle(1:2:end)*pi/180;
    la0=reshape(la0,50,50)';
    la1000=look_angle(2:2:end)*pi/180;
    la1000=reshape(la1000,50,50)';
    la0_ps=griddata(gridX,gridY,la0,ij(:,3),ij(:,2),'linear',{'QJ'});
    la1000_ps=griddata(gridX,gridY,la1000,ij(:,3),ij(:,2),'linear',{'QJ'});
    la=la0_ps+(la1000_ps-la0_ps).*hgt/1000;
    lasavename=['la',num2str(psver)];
    save(lasavename,'la');
end

updir=0;
bperpdir=dir('bperp_*.1.in');
if isempty(bperpdir)
    bperpdir=dir('../bperp_*.1.in');
    updir=1;
end
if length(bperpdir)>0
    bperp_mat=zeros(n_ps,n_image,'single');
    for i=setdiff([1:n_image],master_ix);
        bperp_fname=['bperp_',datestr(day(i),'yyyymmdd'),'.1.in'];
        if updir==1
            bperp_fname=['../',bperp_fname];
        end
        bperp_grid=load(bperp_fname);
        bp0=bperp_grid(1:2:end);
        bp0=reshape(bp0,50,50)';
        bp1000=bperp_grid(2:2:end);
        bp1000=reshape(bp1000,50,50)';
        bp0_ps=griddata(gridX,gridY,bp0,ij(:,3),ij(:,2),'linear',{'QJ'});
        %bp1000_ps=griddata(gridX,gridY,bp1000,ij(:,3),ij(:,2),'linear',{'QJ'});
        %bperp_mat(:,i)=bp0_ps+(bp1000_ps-bp0_ps).*hgt/1000;
        bperp_mat(:,i)=bp0_ps;
    end
    %bperp_mat=[bperp_mat(:,1:master_ix-1),zeros(n_ps,1),bperp_mat(:,master_ix:end)]; % insert master-master bperp (zero)
    bperp_mat=bperp_mat(:,ifgday_ix(:,2))-bperp_mat(:,ifgday_ix(:,1));
else
    bperp_mat=repmat(single(bperp)',n_ps,1);
end

bpsavename=['bp',num2str(psver)];
save(bpsavename,'bperp_mat');

end % end-if



