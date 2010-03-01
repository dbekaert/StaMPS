function [demsavename]=ps_load_dem()
%PS_LOAD_DEM Load DEM, subsample, and save for plotting purposes
%
%   Andy Hooper, July 2006
%  
%   ======================================================
%   08/2006 AH: amended for patches
%   ======================================================


min_posting=0.0001388888889*3; % minimum posting of dem in degrees

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
  dem_lon=str2num(fgetl(fid)); % upper left corner
  dem_lat=str2num(fgetl(fid)); % ditto
  dem_posting=str2num(fgetl(fid));
  dem_format=fgetl(fid);
  fclose(fid);
 
  dem_lat=dem_lat - (dem_length-1)*dem_posting; % lower left corner

  fid=fopen(demname,'r');
  if strfind(upper(dem_format),'I2')
      dem=fread(fid,[dem_width,inf],'int16');
  else    
      dem=fread(fid,[dem_width,inf],'float');
  end    
  dem=dem';
  fclose(fid);
  posting=llh2local(ll0'+[dem_posting; dem_posting]*(dem_width/4),ll0)/(dem_width/4)*1000;
  x_posting=posting(1);
  y_posting=posting(2);
  x_offset=(ll0(1)-dem_lon)/dem_posting*x_posting;
  y_offset=(ll0(2)-dem_lat)/dem_posting*y_posting;

  if dem_posting<min_posting
    scale_factor=1/(ceil(min_posting/dem_posting));
    dem=imresize(dem,scale_factor,'bicubic');
    x_posting=x_posting/scale_factor;
    y_posting=y_posting/scale_factor;
    [dem_length,dem_width]=size(dem);
    dem_lon=dem_lon+(1/scale_factor-1)/2*dem_posting;
    dem_lat=dem_lat+(1/scale_factor-1)/2*dem_posting;
    dem_posting=dem_posting/scale_factor;
  end    

  save(demsavename,'dem','x_offset', 'y_offset', 'x_posting', 'y_posting', 'dem_lon', 'dem_lat', 'dem_width', 'dem_length', 'dem_posting')

%end % end-if



