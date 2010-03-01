function [dd,ar]=combine_amp_dem(az_down,rg_right,exponent,scale_factor,use_slope)
%COMBINE_AMP_DEM combine amp and dem and display
%    MAKE_AMP_DEM(BLUE_DOWN,BLUE_RIGHT,RED_CONTRAST,BLUE_BRIGHTNESS)
%    RED_CONTRAST  (default 0.5) = contrast of interferogram amplitude
%    BLUE_BRIGHTNESS (default 1) = brightness of DEM or DEM slope
%    USE_SLOPE (default 1) = set to 0 to plot DEM elevation instead 
%
%   Modified from make_amp_dem, Andy Hooper, Mar 2005
%
%   ===========================================================
%   05/2007 AH: memory needs reduced  
%   10/2007 AH: added DEM slope option (default)
%   03/2008 AH: memory needs further reduced
%   05/2008 AH: check if multilooked amp file needs updating
%   08/2008 AH: calculate dem_radar.raw format from file length
%   10/2008 AH: Generalise for other sensors
%   10/2009 AH: memory needs reduced  
%   ===========================================================


if nargin<1
    az_down=0;
end

if nargin<2
    rg_right=0;
end

if nargin<3
    exponent=0.5;
end

if nargin<4
    scale_factor=1;
end

if nargin<5
    use_slope=1;
end


%!grep Data_output_format dem.out | gawk 'END {print $2}' > demformat.txt
%fid=fopen('demformat.txt');
%dem_format=fgetl(fid);
%fclose(fid);

!grep 'Number of lines (multilooked)' dem.out | gawk 'END {print $5}' > ml_length.txt
ny=load('ml_length.txt');
!grep 'Number of pixels (multilooked)' dem.out | gawk 'END {print $5}' > width.txt
nx=load('width.txt');
!grep 'Multilookfactor_azimuth_direction' dem.out | gawk 'END {print $2}' > mlaz.txt
ar=load('mlaz.txt');

fid=fopen('dem_radar.raw'); 
fseek(fid,0,1);
fbytes=ftell(fid);
fseek(fid,0,-1);
pbytes=fbytes/ny/nx;
if pbytes==4
    fprintf('dem_radar.raw looks to be single precision\n')
    bb=fread(fid,[nx,inf],'float=>single');
elseif pbytes==8
    fprintf('dem_radar.raw looks to be double precision\n')
    bb=fread(fid,[nx,inf],'double=>single');
else
    fprintf('dem_radar.raw length is %d\n',fbytes)
    fprintf('Number of lines (multilooked) is %d\n',ny)
    fprintf('Number of pixels (multilooked) is %d\n',nx)
    error('dem_radar.raw length is apparently incorrect')
end

%if strfind(upper(dem_format),'I2')
%    bb=fread(fid,[width,inf],'int16=>single');
%else
%    bb=fread(fid,[width,inf],'float=>single');
%end
fclose(fid);
bb=bb';
if use_slope~=0
    bb=[zeros(size(bb,1),1),diff(bb,1,2)];
    scale_factor=scale_factor*100;
end
bb(bb<0)=0;
%[ny,nx]=size(bb);

ampfile=['amp_',num2str(ar),'laz.flt'];

if exist(ampfile,'file')
    ampdir=dir(ampfile);
    cintdir=dir('cint.raw');
    if isfield(ampdir,'datenum')
        ampdatenum=ampdir.datenum;
        cintdatenum=cintdir.datenum;
    else
        ampdatenum=datenum(ampdir.date);
        cintdatenum=datenum(cintdir.date);
    end
    if cintdatenum < ampdatenum
        fid=fopen(ampfile); 
        aa=fread(fid,[nx,inf],'float=>single');
    end
end

if ~exist('aa','var')
    fid=fopen('cint.raw'); 
    aa=zeros(nx,ny,'single'); 
    for i=1:ny
        [aa_bit,byte_count]=fread(fid,[nx*2,ar],'float=>single');
        aa(:,i)=sum(abs(complex(aa_bit(1:2:end,:),aa_bit(2:2:end,:))),2);
    end
    %aa=fread(fid,[width*2,inf],'float=>single');
    %aa=complex(aa(1:2:end,:),aa(2:2:end,:));
    %aa=abs(aa);
    %if size(aa,2)~=size(bb,1)

    %    aa=aa(:,1:5:end)+aa(:,2:5:end)+aa(:,3:5:end)+aa(:,4:5:end)+aa(:,5:5:end);
    %end
    fclose(fid);
    fid=fopen(ampfile,'w'); 
    fwrite(fid,aa,'float');
    fclose(fid);
end


cc=zeros(size(bb),'single');
if az_down>0
   i2l=1;
   i2h=ny-az_down;
   i3l=1+az_down;
   i3h=ny;
else
   i2l=1-az_down;
   i2h=ny;
   i3l=1;
   i3h=ny+az_down;
end

if rg_right>0
   j2l=1;
   j2h=nx-rg_right;
   j3l=1+rg_right;
   j3h=nx;
else
   j2l=1-rg_right;
   j2h=nx;
   j3l=1;
   j3h=nx+rg_right;
end

cc(i3l:i3h,j3l:j3h)=bb(i2l:i2h,j2l:j2h);
clear bb
%cc_sort=sort(cc(:));
%cc(cc<cc_sort(round(nx*ny*0.5)))=cc(round(cc_sort(round(nx*ny*0.5))));
%hist(cc(:),100)
aa=aa';
%aa=log(aa);
aa=aa.^exponent;
aa_sort=sort(aa(:));
aa_max=aa_sort(round(length(aa_sort)*0.99));
aa=aa-min(aa(:));
aa=aa./aa_max*255;
cc=cc-min(cc(:));
cc=cc./max(cc(:))*255*scale_factor;
dd=zeros(size(cc,1),size(cc,2),3,'uint8');
dd(:,:,1)=aa;
dd(:,:,3)=cc;
