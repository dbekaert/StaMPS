function [demsavename]=ps_load_dem(lon1,lon2,lat1,lat2)
%PS_LOAD_DEM Load DEM, subsample, and save for plotting purposes
%   PS_LOAD_DEM(LON1,LON2,LAT1,LAT2) 
%   where inputs specify optional crop area
%
%   Andy Hooper, July 2006
%  
%   ======================================================
%   08/2006 AH: amended for patches
%   02/2010 AH: add option to set posting
%   09/2010 AH: make pixels square
%   09/2010 AH: crop to area of selected pixels
%   06/2016 AH: allow for different posting in x and y 
%   ======================================================

if nargin<4
    load psver
    psname=['ps',num2str(psver)];
    ps=load(psname);
    ll_border=((max(ps.lonlat)-min(ps.lonlat))*0.05);
end

if nargin<1
   lon1=min(ps.lonlat(:,1))-ll_border(1);
end
if nargin<2
   lon2=max(ps.lonlat(:,1))+ll_border(1);
end
if nargin<3
   lat1=min(ps.lonlat(:,2))-ll_border(2);
end
if nargin<4
   lat2=max(ps.lonlat(:,2))+ll_border(2);
end



%min_posting=getparm('plot_dem_posting',1)*9e-6; %minimum posting of dem in degrees
min_posting=getparm('plot_dem_posting',1); %minimum posting of dem 

load psver
psname=['ps',num2str(psver)];
ps=load(psname);
if isfield(ps,'ll0')
    ll0=ps.ll0;
else
    ll0=mean(ps.lonlat);
end
clear ps

demparmsname='./demparms.in';
demsavename=['./dem'];
if ~exist(demparmsname,'file')
    demparmsname='../demparms.in';
    demsavename=['../dem'];
end
if ~exist(demparmsname,'file')
    demparmsname='../../demparms.in';
    demsavename=['../../dem'];
end
if ~exist(demparmsname,'file')
    error('Cannot find demparms.in in this directory, parent or grandparent')
end

fid=fopen(demparmsname,'r');
%if fid>0
  demname=fgetl(fid);
  dem_width=str2num(fgetl(fid));
  dem_length=str2num(fgetl(fid));
  dem_lon1=str2num(fgetl(fid)); % upper left corner
  dem_lat2=str2num(fgetl(fid)); % ditto
  dem_posting=str2num(fgetl(fid)); % can be 1 or 2 numbers
  dem_posting=fliplr(dem_posting); % put lon posting first (if there are 2)
  dem_format=fgetl(fid);
  fclose(fid);
 
  dem_lat1=dem_lat2 - (dem_length-1)*dem_posting(end); % lower right corner
  dem_lon2=dem_lon1 + (dem_width-1)*dem_posting(1); % lower right corner

  fid=fopen(demname,'r');
  if strfind(upper(dem_format),'I2')
      dem=fread(fid,[dem_width,inf],'int16');
  else    
      dem=fread(fid,[dem_width,inf],'float');
  end    
  dem=dem';
  fclose(fid);

  i1=1;
  i2=dem_length;
  j1=1;
  j2=dem_width;

  if dem_lon1<lon1
      j1=round((lon1-dem_lon1)/dem_posting(1))+1;
  end
  if dem_lon2>lon2
      j2=round((lon2-dem_lon1)/dem_posting(1))+1;
  end
  if dem_lat1<lat1
      i2=round((dem_lat2-lat1)/dem_posting(end))+1;
  end
  if dem_lat2>lat2
      i1=round((dem_lat2-lat2)/dem_posting(end))+1;
  end

  dem=dem(i1:i2,j1:j2);
  dem_lon1=dem_lon1+(j1-1)*dem_posting(1);
  dem_lat1=dem_lat2-(i2-1)*dem_posting(end);
  dem_length=i2-i1+1;
  dem_width=j2-j1+1;

  posting=llh2local(ll0'+[dem_posting(1); dem_posting(end)]*(dem_width/4),ll0)/(dem_width/4)*1000;
  x_posting=posting(1);
  y_posting=posting(2);
  x_offset=(ll0(1)-dem_lon1)/dem_posting(1)*x_posting;
  y_offset=(ll0(2)-dem_lat1)/dem_posting(end)*y_posting;

  if max(posting)<min_posting
    %scale_factor=1/(ceil(min_posting/dem_posting));
    i_new=round(dem_length*y_posting/min_posting);
    j_new=round(dem_width*x_posting/min_posting);
    %dem=imresize(dem,scale_factor,'bicubic');
    dem=imresize(dem,[i_new,j_new],'bicubic');
    scale_factor_x=dem_width/j_new;
    scale_factor_y=dem_length/i_new;
    x_posting=x_posting*scale_factor_x;
    y_posting=y_posting*scale_factor_y;
    [dem_length,dem_width]=size(dem);
    dem_lon=dem_lon1+(1/scale_factor_x-1)/2*dem_posting(1);
    dem_lat=dem_lat1+(1/scale_factor_y-1)/2*dem_posting(end);
    dem_posting_lon=dem_posting(1)*scale_factor_x;
    dem_posting_lat=dem_posting(end)*scale_factor_y;
  else
    dem_lon=dem_lon1;
    dem_lat=dem_lat1;
    dem_posting_lon=dem_posting(1);
    dem_posting_lat=dem_posting(2);
    setparm('plot_dem_posting',round(max(posting)))
  end    

  save(demsavename,'dem','x_offset', 'y_offset', 'x_posting', 'y_posting', 'dem_lon', 'dem_lat', 'dem_width', 'dem_length', 'dem_posting_lon','dem_posting_lat');

%end % end-if



