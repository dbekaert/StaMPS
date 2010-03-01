function [ph_lims]=ps_plot_ifg(in_ph,bg_flag,col_rg,lon_rg,lat_rg)
% PS_PLOT_IFG plot phase of a PS interferogram
%    PS_PLOT_IFG(PHASE,BACKGROUND,LIMS)
%
%    PHASE is a vector of phases (wrapped or unwrapped)
%    
%    BG_FLAG,  0, black background, lon/lat axes (default)
%              1, white background, lon/lat axes
%              2, shaded relief topo, lon/lat axes
%              3, 3D topo, lon/lat axes
%              4, mean amplitude image
%              5, mean amplitude image, brightness showing through PS
%              6, white background, xy axis (rotated lon/lat)
%
%    COL_RG is the color map limits (defaults to use maximum range)  
%    LON_RG is the longitude range (defaults to use maximum range)            
%    LAT_RG is the latitude range (defaults to use maximum range)            
%
%   Andy Hooper, June 2006

%   ======================================================================
%   01/2010 AH: Fix bug for amplitude plots in PATCH directories
%   11/2009 AH: Fix bug for BG_FLAG 6 to work
%   11/2009 AH: Reduce memory needs for amplitude plots
%   11/2009 AH: Fix orientation and aspect ratio for amplitude plots
%   ======================================================================



plot_pixel_size=getparm('plot_pixel_size');
%pixel_aspect_ratio=getparm('pixel_aspect_ratio');
plot_color_scheme=getparm('plot_color_scheme');
shade_rel_angle=getparm('shade_rel_angle');
lonlat_offset=getparm('lonlat_offset');
heading=getparm('heading');

if nargin < 1
    error('PS_PLOT_IFG(PHASE,BACKGROUND,LIMS)')
end

if nargin < 3
    col_rg=[];
end

if nargin < 4
    lon_rg=[];
end

if nargin < 5
    lat_rg=[];
end

if nargin < 2
    bg_flag=1; % background flag
end

if bg_flag==4 | bg_flag==5
    bg_amp=1;
else
    bg_amp=0;
end

marker_size=8; 


x_posting=4/3600;
y_posting=4/3600;

load psver
psname=['ps',num2str(psver)];

ps=load(psname);
lonlat=ps.lonlat;
if bg_amp==1
    ij=ps.ij;
end

lonlat(:,1)=lonlat(:,1)+lonlat_offset(1);
lonlat(:,2)=lonlat(:,2)+lonlat_offset(2);

if ~isempty(lon_rg)
    ix=lonlat(:,1)>=lon_rg(1)&lonlat(:,1)<=lon_rg(2);
    in_ph=in_ph(ix);
    lonlat=lonlat(ix,:);
    if bg_amp==1
        ij=ij(ix,:);
    end
end

if ~isempty(lat_rg)
    ix=lonlat(:,2)>=lat_rg(1)&lonlat(:,2)<=lat_rg(2);
    in_ph=in_ph(ix);
    lonlat=lonlat(ix,:);
    if bg_amp==1
        ij=ij(ix,:);
    end
end

min_ph=0;
max_ph=0;

if ~isempty(in_ph)
if ~isempty(col_rg)
    min_ph=min(col_rg);
    max_ph=max(col_rg);
else
    if isreal(in_ph)
        min_ph=min(in_ph);
        max_ph=max(in_ph);
    end
end
ph_range=max_ph-min_ph;

if ~isreal(in_ph) | ph_range==0
    if ph_range==0
        min_ph=-pi;
        max_ph=+pi;
        ph_range=2*pi;
    end
    in_ph=angle(in_ph);
    c=hsv(64);
else
    if strncmpi(plot_color_scheme,'inflation',9)
        c=flipud(jet(64));
    else
        c=jet(64);
    end
end

col_ix=round(((in_ph-min_ph)*63/ph_range)+1);

col_ix(col_ix>64)=64;
col_ix(col_ix<1)=1;
else
    c=hsv(64);
end
hold on

if bg_flag==4 % plot on amplitude image

    ampfile='./amp_mean.mat';
    if ~exist(ampfile,'file')
        ampile='../amp_mean.mat';
    end
    if ~exist(ampfile,'file')
        ampfile=ps_load_mean_amp;
    end
    load(ampfile)

    if exist('./patch.in','file')
       patch=load('./patch.in');
       start_rg=patch(1);
       start_az=patch(3);
       ij(:,2)=ij(:,2)-start_az+1;
       ij(:,3)=ij(:,3)-start_rg+1;
    end
    
    amp_mean=uint16(amp_mean);
    [n,m]=size(amp_mean);
    c=[gray(256);c];

    mean_x_post=mean((ps.xy(:,2)-mean(ps.xy(:,2)))./(ps.ij(:,3)-mean(ps.ij(:,3))));
    mean_y_post=mean((ps.xy(:,3)-mean(ps.xy(:,3)))./(ps.ij(:,2)-mean(ps.ij(:,2))));
    pixel_aspect_ratio=abs(mean_x_post/mean_y_post);

    az_ix=pixel_aspect_ratio:pixel_aspect_ratio:m*pixel_aspect_ratio;
    
    pixel_margin1=floor((plot_pixel_size-1)/2);
    pixel_margin2=ceil((plot_pixel_size-1)/2);

    for i=1 : length(in_ph)
        ix1=ij(i,2)-floor((pixel_aspect_ratio*plot_pixel_size/2)-1):ij(i,2)+ceil(pixel_aspect_ratio*plot_pixel_size/2);
        ix1=ix1(ix1>0&ix1<=n);
        ix2=ij(i,3)-pixel_margin1+1:ij(i,3)+pixel_margin2+1;
        ix2=ix2(ix2>0&ix2<=m);
        if ~(isnan(col_ix(i)))
            amp_mean(ix1,ix2)=256+col_ix(i);
        end
    end
   
    if cos(heading*pi/180)<0 
        amp_mean=flipud(fliplr(amp_mean));
    end 

    image(az_ix,[1:size(amp_mean,1)],amp_mean)
    set(gca,'xtick',[],'ytick',[])
    
    axis equal
    axis tight
    
elseif bg_flag==5 % plot on amplitude image, let amp show through color
    
    ampfile='./amp_mean.mat';
    if ~exist(ampfile,'file')
        ampfile='../amp_mean.mat';
    end
    if ~exist(ampfile,'file')
        ampfile=ps_load_mean_amp;
    end
    load(ampfile)

    if exist('./patch.in','file')
       patch=load('./patch.in');
       start_rg=patch(1);
       start_az=patch(3);
       ij(:,2)=ij(:,2)-start_az+1;
       ij(:,3)=ij(:,3)-start_rg+1;
    end

    n_gray=16;
    n_color=64;
    
    ci=single(amp_mean);
    clear amp_mean
    ci=ci-min(ci(:));
    cimax=max(ci(:))*(1+1e-6);
    ci=uint16(floor(ci/cimax*n_gray)*(n_color+1));
    c_base=c;
    c=zeros(n_gray*(n_color+1),3);
    for i=1:n_gray
        i_frac=i/n_gray;
        c((i-1)*(n_color+1)+1:i*(n_color+1),:)=[i_frac i_frac i_frac;c_base*(0.6+i_frac*0.4)];
    end
    cd=ci;
    [n,m]=size(cd);

    mean_x_post=mean((ps.xy(:,2)-mean(ps.xy(:,2)))./(ps.ij(:,3)-mean(ps.ij(:,3))));
    mean_y_post=mean((ps.xy(:,3)-mean(ps.xy(:,3)))./(ps.ij(:,2)-mean(ps.ij(:,2))));
    pixel_aspect_ratio=abs(mean_x_post/mean_y_post);
    az_ix=pixel_aspect_ratio:pixel_aspect_ratio:m*pixel_aspect_ratio;
    
    pixel_margin1=floor((plot_pixel_size-1)/2);
    pixel_margin2=ceil((plot_pixel_size-1)/2);

    for i=1 : length(in_ph)
        ix1=ij(i,2)-floor((pixel_aspect_ratio*plot_pixel_size/2)-1):ij(i,2)+ceil(pixel_aspect_ratio*plot_pixel_size/2);
        ix1=ix1(ix1>0&ix1<=n);
        ix2=ij(i,3)-pixel_margin1+1:ij(i,3)+pixel_margin2+1;
        ix2=ix2(ix2>0&ix2<=m);
        if ~(isnan(col_ix(i)))
            cd(ix1,ix2)=ci(ix1,ix2)+uint16(col_ix(i));
        end
    end
   
    %if abs(heading)>90 
    %    amp_mean=flipud(fliplr(cd));
    %end 
    if cos(heading*pi/180)<0 
        amp_mean=flipud(fliplr(cd));
    end 

    image(az_ix,[1:size(cd,1)],cd)
    set(gca,'xtick',[],'ytick',[])
    
    axis equal
    axis tight
    
elseif floor(bg_flag)==2
    
    demfile='dem.mat';
    if ~exist(demfile,'file')
        demfile='../dem.mat';
    end
    if ~exist(demfile,'file')
        demfile='../../dem.mat';
    end
    if ~exist(demfile,'file')
        demfile=ps_load_dem;
    end
    load(demfile)

    [dem_y,dem_x]=size(dem);
    c2=gray(64);
    c2=c2(35:50,:);
    x=[dem_lon:dem_posting:(dem_x-1)*dem_posting+dem_lon];
    y=[(dem_y-1)*dem_posting+dem_lat:-dem_posting:dem_lat];
    demx=round((lonlat(:,1)-dem_lon)/dem_posting)+1;
    demy=round((y(1)-lonlat(:,2))/dem_posting)+1;
    [X,Y]=meshgrid(x,y);
    [cindx,cimap,clim] = shaderel(X,Y,dem,[0.7 0.7 0.7],shade_rel_angle);
    if bg_flag==2
        c=[cimap;0,0,0;c];
    else
        c=[cimap;1,1,1;c];
    end
    cindx(dem==0)=257;
    h=image(x,y,cindx);
    pixel_margin1=floor((plot_pixel_size-1)/2);
    pixel_margin2=ceil((plot_pixel_size-1)/2);
    R=get(h,'cdata');
    pixel_margin1=floor((plot_pixel_size-1)/2);
    pixel_margin2=ceil((plot_pixel_size-1)/2);
    
   	for i=1 : length(in_ph)
        %if ~(isnan(col_ix(i))) & demy(i)>0 & demy(i)<=size(X,1) & demx(i)>0 & demx(i)<=size(X,2) & dem(demy(i),demx(i))~=0
        if ~(isnan(col_ix(i))) & demy(i)>0 & demy(i)<=size(X,1) & demx(i)>0 & demx(i)<=size(X,2) 
            ix1=demy(i)-pixel_margin1:demy(i)+pixel_margin2;
            ix2=demx(i)-pixel_margin1:demx(i)+pixel_margin2;
            ix1=ix1(ix1>0&ix1<=size(X,1));
            ix2=ix2(ix2>0&ix2<=size(X,2));
            R(ix1,ix2)=col_ix(i)+257;
        end
    end

    set(h,'cdata',R);

    axis tight
    dem_length=dem_y;
    dem_width=dem_x;
    lat_range=dem_posting*dem_length;
    lon_range=dem_posting*dem_width;
    xy_ratio=llh2local([dem_lon+lon_range;dem_lat+lat_range;0],[dem_lon;dem_lat;0]);
    aspect_ratio=[xy_ratio(1)/xy_ratio(2),1,1];
    set(gca,'plotboxaspectratio',aspect_ratio)

    
elseif bg_flag==3    % plot on 3D DEM
    
    
    demfile='dem.mat';
    if ~exist(demfile,'file')
        demfile='../dem.mat';
    end
    if ~exist(demfile,'file')
        demfile='../../dem.mat';
    end
    if ~exist(demfile,'file')
        demfile=ps_load_dem;
    end
    load(demfile)
    
    [dem_y,dem_x]=size(dem);
    c2=gray(80);
    x=[dem_lon:dem_posting:(dem_x-1)*dem_posting+dem_lon];
    y=[(dem_y-1)*dem_posting+dem_lat:-dem_posting:dem_lat];

    [X,Y]=meshgrid(x,y);
    c=[c2(8:71,:);c];
    
    demx=round((lonlat(:,1)-dem_lon)/dem_posting)+1;
    demy=round((y(1)-lonlat(:,2))/dem_posting)+1;
 
    h=surfl(X,Y,dem,shade_rel_angle);
    dem_length=dem_y;
    dem_width=dem_x;
    lat_range=dem_posting*dem_length;
    lon_range=dem_posting*dem_width;
    xy_ratio=llh2local([dem_lon+lon_range;dem_lat+lat_range;0],[dem_lon;dem_lat;0]);
    aspect_ratio=[xy_ratio(1)/xy_ratio(2),1,0.1];
    set(gca,'plotboxaspectratio',aspect_ratio)
    
    shading interp
    R=get(h,'cdata');
    R=R/2;
    
    c=[c;0,0,0];
    R(dem==0)=129;

    pixel_margin1=floor((plot_pixel_size-1)/2);
    pixel_margin2=ceil((plot_pixel_size-1)/2);
    
   	for i=1 : length(in_ph)
        if ~(isnan(col_ix(i)))
            ix1=demy(i)-pixel_margin1:demy(i)+pixel_margin2;
            ix2=demx(i)-pixel_margin1:demx(i)+pixel_margin2;
            ix1=ix1(ix1>0&ix1<=size(X,1));
            ix2=ix2(ix2>0&ix2<=size(X,2));
            R(ix1,ix2)=col_ix(i)/128+0.5;
        end
    end

    set(h,'cdata',R);
    view(20,40)
    axis off
    
    
elseif bg_flag==0 | bg_flag==1    % lon/lat axes
    
    
    x=[min(lonlat(:,1))-2*x_posting:x_posting:max(lonlat(:,1))+2*x_posting];
    y=[min(lonlat(:,2))-2*y_posting:y_posting:max(lonlat(:,2))+2*y_posting];
    
    [X,Y]=meshgrid(x,y);
    if bg_flag==0
        c=[[0 0 0];c]; % black background
    else
        c=[[1 1 1];c]; % white background
    end
    
    demx=round((lonlat(:,1)-x(1))/x_posting)+1;
    demy=round((lonlat(:,2)-y(1))/y_posting)+1;

    R=zeros(size(X));
    
    pixel_margin1=floor((plot_pixel_size-1)/2);
    pixel_margin2=ceil((plot_pixel_size-1)/2);
    
   	for i=1 : length(in_ph)
        if ~(isnan(col_ix(i)))
            ix1=demy(i)-pixel_margin1:demy(i)+pixel_margin2;
            ix2=demx(i)-pixel_margin1:demx(i)+pixel_margin2;
            ix1=ix1(ix1>0&ix1<=size(X,1));
            ix2=ix2(ix2>0&ix2<=size(X,2));
            R(ix1,ix2)=col_ix(i)+1;
        end
    end
    
    dem_length=size(R,1);
    dem_width=size(R,2);
    lat_range=y_posting*dem_length;
    lon_range=x_posting*dem_width;
    dem_lat=y(1);
    dem_lon=x(1);
    xy_ratio=llh2local([dem_lon+lon_range;dem_lat+lat_range;0],[dem_lon;dem_lat;0]);
    aspect_ratio=[xy_ratio(1)/xy_ratio(2),1,1];
    set(gca,'plotboxaspectratio',aspect_ratio)
    
    image(x,y,R)
    axis tight
    
   
    

elseif bg_flag==6     % xy axes
    
    x_posting=10;
    y_posting=10;
    
    x=[min(ps.xy(:,2))-2*x_posting:x_posting:max(ps.xy(:,2))+2*x_posting];
    y=[min(ps.xy(:,3))-2*y_posting:y_posting:max(ps.xy(:,3))+2*y_posting];
    
    [X,Y]=meshgrid(x,y);
    c=[[1 1 1];c]; % white background
    
    demx=round((ps.xy(:,2)-x(1))/x_posting)+1;
    demy=round((ps.xy(:,3)-y(1))/y_posting)+1;

    R=zeros(size(X));
    
    pixel_margin1=floor((plot_pixel_size-1)/2);
    pixel_margin2=ceil((plot_pixel_size-1)/2);
    
   	for i=1 : length(in_ph)
        if ~(isnan(col_ix(i)))
            ix1=demy(i)-pixel_margin1:demy(i)+pixel_margin2;
            ix2=demx(i)-pixel_margin1:demx(i)+pixel_margin2;
            ix1=ix1(ix1>0&ix1<=size(X,1));
            ix2=ix2(ix2>0&ix2<=size(X,2));
            R(ix1,ix2)=col_ix(i)+1;
        end
    end
    
    image(x,y,R)
    set(gca,'xtick',[],'ytick',[])
    axis ij
    axis equal
    axis tight
    
    

else    
    
 	for i=1 : length(in_ph)
        if ~(isnan(col_ix(i)))
            p=plot(ps.xy(i,2),ps.xy(i,3),'.');
		    set(p,'color',c(col_ix(i),:));   
        end
    end
    axis equal
    axis tight
end

ph_lims=[max_ph,min_ph];

hold off
colormap(c);



