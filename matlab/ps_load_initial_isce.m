function []=ps_load_initial_isce(data_inc)
%PS_LOAD_INITIAL Initial load of PS files into matlab workspace
%   PS_LOAD_INITIAL() loads various files and stores them
%   in various workspaces. The version number, PSVER, is set to 1.
%
%   Andy Hooper, June 2006
%  
%   ======================================================================
%   06/2006 AH: Initialize day variable
%   07/2006 AH: Fix error in la and bperp matrix calculation
%   07/2006 AH: Add patch compatibility
%   09/2006 AH: Store wrapped phase in separate workspace
%   11/2007 AH: Change bperp to be the bperp at 0 m
%   01/2008 AH: Selection of bperp files made more rigourous
%   03/2009 AH: Change bperp back to be the bperp at local height
%   10/2009 AH: Add n_image to mat file
%   04/2010 KS: Fixed small issue with new logit function
%   10/2009 MA: Oversampling factor file introduced
%   12/2012 DB: Add compatiblility with Matlab2012B, keep backward compatible
%   01/2016 DB: Replace save with stamps_save which checks for var size when saving 
%   04/2017 DB: Include fix for isce2stamps calamp.out slc filenaming
%   04/2018 DB: Add support for baseline grids
%   ======================================================================


%NB IFGS assumed in ascending date order
fprintf('Loading data into matlab...\n')

phname=['pscands.1.ph'];            % for each PS candidate, a float complex value for each ifg
ijname=['pscands.1.ij'];            % ID# Azimuth# Range# 1 line per PS candidate
bperpname=['bperp.1.in'];           % in meters 1 line per ifg
dayname=['day.1.in'];               % YYYYMMDD, 1 line per ifg
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
incname=['inc_angle.raw'];
incsavename=['inc',num2str(psver)];
lasavename=['la',num2str(psver)];

% Get matlab version as function arguments change with the matlab version
matlab_version = version('-release');           % [DB] getting the matlab version
matlab_version = str2num(matlab_version(1:4));  % [DB] the year

if ~exist(dayname,'file')
    dayname= ['../',dayname];
end
day=load(dayname);
year=floor(day/10000);
month=floor((day-year*10000)/100);
monthday=day-year*10000-month*100;
slave_day=datenum(year,month,monthday);
[slave_day,day_ix]=sort(slave_day);

if ~exist(masterdayname,'file')
    masterdayname= ['../',masterdayname];
end
master_day=load(masterdayname);
master_day_yyyymmdd=master_day;
year=floor(master_day/10000);
month=floor((master_day-year*10000)/100);
monthday=master_day-year*10000-month*100;
master_day=datenum(year,month,monthday);

master_ix=sum(slave_day<master_day)+1;
day=[slave_day(1:master_ix-1);master_day;slave_day(master_ix:end)]; % insert master day


%% bperp one value per slave
if ~exist(bperpname,'file')
    bperpname= ['../',bperpname];
end
bperp=load(bperpname); 
bperp=bperp(day_ix);
bperp=[bperp(1:master_ix-1);0;bperp(master_ix:end)]; % insert master-master bperp (zero)
n_ifg=size(bperp,1);
n_image=n_ifg;

%% heading 
if ~exist(headingname,'file')
    headingname= ['../',headingname];
end
heading=load(headingname);
if isempty(heading)
    error('heading.1.in is empty')
end
setparm('heading',heading,1);

%% wavelength
if ~exist(lambdaname,'file')
    lambdaname= ['../',lambdaname];
end
lambda=load(lambdaname);
setparm('lambda',lambda,1);

%% Radar coordinates
ij=load(ijname);
n_ps=size(ij,1);

if ~exist(calname,'file')
    calname= ['../',calname];
end
if exist(calname,'file')
    [calfile,calconst]=textread(calname,'%s%f');
    caldate=zeros(length(calfile),1);
    for i = 1 : length(calfile)
        aa=strread(calfile{i},'%s','delimiter','/');
        try 
            bb=str2num(aa{end}(1:8));
            if isempty(bb)
                bb=str2num(aa{end-1}(1:8));
            end
        catch
            if strcmpi(aa{end-1},'master');
                bb=str2num(aa{end-2}(end-7:end));
            end            
        end
        caldate(i)=bb;    
    end
    not_master_ix=caldate~=master_day_yyyymmdd;
    caldate=caldate(not_master_ix);
    calconst=calconst(not_master_ix);
    [Y,I]=sort(caldate);
    calconst=calconst(I);
else
    calconst=ones(n_ifg-1,1);
end

fid=fopen(phname,'r');
ph=zeros(n_ps,n_ifg-1,'single');
byte_count=n_ps*2;
for i=1:n_ifg-1
    [ph_bit,byte_count]=fread(fid,[(n_ps)*2,1],'float');
    ph_bit=single(ph_bit);
    ph(:,i)=complex(ph_bit(1:2:end),ph_bit(2:2:end));
end

ph=ph(:,day_ix);
zero_ph=sum(ph==0,2);
nonzero_ix=zero_ph<=1;       % if more than 1 phase is zero, drop node
ph=ph./repmat(calconst',n_ps,1);  % scale amplitudes
ph=[ph(:,1:master_ix-1),ones(n_ps,1),ph(:,master_ix:end)]; % insert zero phase master-master ifg


%% LON LAT
fid=fopen(llname,'r');
lonlat=fread(fid,[2,inf],'float');
lonlat=lonlat';
fclose(fid);

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
stamps_save(savename,ij,lonlat,xy,bperp,day,master_day,master_ix,n_ifg,n_image,n_ps,sort_ix,ll0,calconst,master_ix,day_ix);
save psver psver

phsavename=['ph',num2str(psver)];
% save(phsavename,'ph');
stamps_save(phsavename,ph);

if exist(daname,'file')
  D_A=load(daname);
  D_A=D_A(sort_ix);
  dasavename=['da',num2str(psver)];
  stamps_save(dasavename,D_A);
end


if exist(hgtname,'file')
    fid=fopen(hgtname,'r');
    hgt=fread(fid,[1,inf],'float');
    hgt=hgt';
    hgt=hgt(sort_ix);
    fclose(fid);
    hgtsavename=['hgt',num2str(psver)];
    stamps_save(hgtsavename,hgt);
end

if ~exist(widthname,'file')
    widthname= ['../',widthname];
end
width=load(widthname);

if ~exist(lenname,'file')
    lenname= ['../',lenname];
end
len=load(lenname);



% try to see if inc angle exist, fall back to la if needed
if nargin<1 | isempty(data_inc)
    clear data_inc
    % might be in dirs above when called from patch dir
    % inc angle
    if ~exist(incname,'file')
        incname= ['..',filesep,incname];
        if ~exist(incname,'file')
            incname= ['..',filesep,incname];
        end
    end
    
    if exist(incname,'file')
         data_inc = single(load_isce(incname));
    end
end
% checking if the inc angle exist if not try look angle
if exist('data_inc','var')
    % inc angle exist
    % getting the position of the PS
    IND = sub2ind(size(data_inc),ij(:,3)+1,ij(:,2)+1);
    % storing the data
    inc=data_inc(IND);
    inc = inc*pi./180;
    stamps_save(incsavename,inc)
else
    % trying look angle instead 
     if ~exist(laname,'file')
        laname= ['..',filesep,laname];
        if ~exist(incname,'file')
            laname= ['..',filesep,laname];
            clear laname
        end
     end
     % laod the data if it exists
     if exist(laname,'file')
         data_la = single(load_isce(laname));

        % getting the position of the PS
        IND = sub2ind(size(data_la),ij(:,3)+1,ij(:,2)+1);

        % storing the data
        la=data_la(IND);
        la = la*pi./180;
        stamps_save(lasavename,la)
        clear data_la
     end
end



updir=0;
bperpdir=dir(['..' filesep 'baselineGRID_*']);
if isempty(bperpdir)
    bperpdir=dir(['..' filesep '..' filesep 'baselineGRID_*']);
    updir=1;
end
if length(bperpdir)>0
    % make a set of baseline excluding the master
    bperp_mat=zeros(n_ps,n_image-1,'single');
    counter =1;
    for i=setdiff([1:n_image],master_ix);
        bperp_fname=['..' filesep 'baselineGRID_' datestr(day(i),'yyyymmdd')];
        if updir==1
            bperp_fname=['..' filesep bperp_fname];
        end
        bperp_grid = load_isce(bperp_fname);

        % check if the baseline grid is a full IFG grid or smaller grid
        if size(bperp_grid)==[width, len]
            % getting the position of the PS
            IND = sub2ind(size(bperp_grid),ij(:,3)+1,ij(:,2)+1);
            % storing the data
            bp0_ps=bperp_grid(IND);
            if counter==1
                fprintf('yes \n')
            end
        else
            if counter==1            
                fprintf('no \n')
                [gridX,gridY]=meshgrid(linspace(0,width,size(bperp_grid,1)),linspace(0,len,size(bperp_grid,2)));
            end
            bp0_ps=griddata_version_control(gridX,gridY,bperp_grid',ij(:,3),ij(:,2),'linear',matlab_version);               % [DB] fix matlab2012 version and older
            
        end
        bperp_mat(:,counter)=bp0_ps;
        counter = counter+1;
    end
else
    bperp_mat=repmat(single(bperp([1:master_ix-1,master_ix+1:end]))',n_ps,1);
end
bpsavename=['bp',num2str(psver)];
stamps_save(bpsavename,bperp_mat);

