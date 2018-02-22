function [ph_lims]=ps_plot_ifg(in_ph,bg_flag,col_rg,lon_rg,lat_rg,ext_data)


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
%              7, black background, xy axis (rotated lon/lat)
%
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
%   02/2010 AH: Add options for dem posting and scatterer size
%   02/2010 AH: Plot reference centre if circular
%   09/2010 AH: Add back option to plot scatterer as multiple pixels
%   06/2011 AH: Fix bug so only patch amplitudes plotted for PATCH dirs
%   04/2013 DB: Allow for plotting of additonal data, e.g. LOS GPS
%   05/2013 DB: More flexibility for ext data option
%   05/2013 DB: Fix issues with colorbar for external data plotting option
%   11/2013 DB: Add gray as colorbar option
%   12/2013 DB: Add GMT colorbar options
%   04/2015 DB: Bugfix of the colorbar when memory is not cleared properly.
%               harcoded deflation and inflation using jet, rather than
%               colormap options
%   09/2015 AH: Don't plot dropped patches on sides in amplitude plots
%   12/2015 AH: Speed up plotting when 1 pixel per PS 
%   02/2016 DB: Fix in case the colormatrix has NaN values
%   02/2016 DB: Bug fix for BW background, which made PS pixels as large as the image width. 
%   01/2017 DB: Bug fix in case there is nan for BW background
%   01/2017 DB: Add the hardcoded option to do some filtering
%   03/2017 DB: for lonlat option when you use lon, lat range allow it to fix the axis as well
%   ======================================================================

plot_pixel_m=getparm('plot_scatterer_size');
plot_pixel_size=getparm('plot_pixels_scatterer');
plot_dem_posting=getparm('plot_dem_posting');
plot_color_scheme=getparm('plot_color_scheme');
shade_rel_angle=getparm('shade_rel_angle');
lonlat_offset=getparm('lonlat_offset');
heading=getparm('heading');
ref_radius=getparm('ref_radius');
ref_centre=getparm('ref_centre_lonlat');
small_baseline_flag=getparm('small_baseline_flag');

filter_extra = 'n';
filter_type = 'mean';

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

if nargin<6 || isempty(ext_data)
   plot_ext_data =0;
   ext_data = [];
else
   plot_ext_data =1;
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


load psver
psname=['ps',num2str(psver)];

ps=load(psname);
lonlat=ps.lonlat;
if bg_amp==1
    ij=ps.ij;
end

mean_x_post=(max(ps.xy(:,2))-min(ps.xy(:,2)))/(max(ps.ij(:,3))-min(ps.ij(:,3)));
mean_y_post=(max(ps.xy(:,3))-min(ps.xy(:,3)))/(max(ps.ij(:,2))-min(ps.ij(:,2)));
pixel_aspect_ratio=abs(mean_x_post/mean_y_post);


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
        if strncmpi(plot_color_scheme,'gray',4)
            c=flipud(gray(64));  
        elseif strncmpi(plot_color_scheme,'inflation',9)
            c=flipud(jet(64));
        elseif strncmpi(plot_color_scheme,'GMT_relief',10)
            [c] = cptcmap('GMT_relief','ncol',129) ;    
        	c= c(65:end,:);
        elseif strncmpi(plot_color_scheme,'GMT_globe',9)
            [c] = cptcmap('GMT_globe','ncol',129) ;
                c=[c(66,:) ;c(66:end,:)];
	
        else    % deflation
            c=jet(64);
        end
    end

    col_ix=round(((in_ph-min_ph)*63/ph_range)+1);

    col_ix(col_ix>64)=64;
    col_ix(col_ix<1)=1;
else
    c=hsv(64);
end

% duplicating the color for the plotting of extenal data
c_original = c;
hold on


if bg_flag==4 % plot on amplitude image

    ampfile='./amp_mean.mat';
    if ~exist(ampfile,'file') & strcmpi(small_baseline_flag,'y')
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

    plot_pixel_size=round(abs(plot_pixel_m/mean_x_post));
    
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
   
    ix=sum(amp_mean)~=0; % dropped patches
    amp_mean=amp_mean(:,ix);
    m=sum(ix);
    ix=sum(amp_mean,2)~=0; % dropped patches
    amp_mean=amp_mean(ix,:);
    
    az_ix=pixel_aspect_ratio:pixel_aspect_ratio:m*pixel_aspect_ratio;
    
    if cos(heading*pi/180)<0 
        amp_mean=flipud(fliplr(amp_mean));
    end 

    image(az_ix,[1:size(amp_mean,1)],amp_mean)
    set(gca,'xtick',[],'ytick',[])
    
    axis equal
    axis tight
    
elseif bg_flag==5 % plot on amplitude image, let amp show through color
    
    ampfile='./amp_mean.mat';
    if ~exist(ampfile,'file') & strcmpi(small_baseline_flag,'y')
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
        i_frac=((i)/(n_gray)).^2;
        c((i-1)*(n_color+1)+1:i*(n_color+1),:)=[i_frac i_frac i_frac;c_base*(0.4+i_frac*0.6)];
    end
    cd=ci;
    [n,m]=size(cd);

    plot_pixel_size=round(abs(plot_pixel_m/mean_x_post));
    
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
    
    ix=sum(cd)~=0; % dropped patches
    cd=cd(:,ix);
    m=sum(ix);
    ix=sum(cd,2)~=0; % dropped patches
    cd=cd(ix,:);
    
    az_ix=pixel_aspect_ratio:pixel_aspect_ratio:m*pixel_aspect_ratio;
   
    if cos(heading*pi/180)<0 
        cd=flipud(fliplr(cd));
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
        load(demfile)
    else
        load(demfile)
        if round(x_posting/plot_dem_posting*10)~=10 % allow 5% error
            demfile=ps_load_dem;
            load(demfile)
        end
    end

    plot_pixel_size=round((plot_pixel_m/(x_posting)-1)/2)*2+1;
    [dem_y,dem_x]=size(dem);
    c2=gray(64);
    c2=c2(35:50,:);
    x=[dem_lon:dem_posting_lon:(dem_x-1)*dem_posting_lon+dem_lon];
    y=[(dem_y-1)*dem_posting_lat+dem_lat:-dem_posting_lat:dem_lat];
    demx=round((lonlat(:,1)-dem_lon)/dem_posting_lon)+1;
    demy=round((y(1)-lonlat(:,2))/dem_posting_lat)+1;
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
    lat_range=dem_posting_lat*dem_length;
    lon_range=dem_posting_lon*dem_width;
    xy_ratio=llh2local([dem_lon+lon_range;dem_lat+lat_range;0],[dem_lon;dem_lat;0]);
    aspect_ratio=[xy_ratio(1)/xy_ratio(2),1,1];
    set(gca,'plotboxaspectratio',aspect_ratio)
    if ref_radius<inf
        p=plot(ref_centre(1),ref_centre(2),'k*');
        set(p,'markersize',10,'linewidth',1)
    end

    
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
        load(demfile)
        if round(x_posting/plot_dem_posting*10)~=10 % allow 5% error
            demfile=ps_load_dem;
            load(demfile)
        end
    end
    load(demfile)
    
    plot_pixel_size=round((plot_pixel_m/(x_posting)-1)/2)*2+1;
    [dem_y,dem_x]=size(dem);
    c2=gray(80);
    x=[dem_lon:dem_posting_lon:(dem_x-1)*dem_posting_lon+dem_lon];
    y=[(dem_y-1)*dem_posting_lat+dem_lat:-dem_posting_lat:dem_lat];

    [X,Y]=meshgrid(x,y);
    c=[c2(8:71,:);c];
    
    demx=round((lonlat(:,1)-dem_lon)/dem_posting_lon)+1;
    demy=round((y(1)-lonlat(:,2))/dem_posting_lat)+1;
 
    h=surfl(X,Y,dem,shade_rel_angle);
    dem_length=dem_y;
    dem_width=dem_x;
    lat_range=dem_posting_lat*dem_length;
    lon_range=dem_posting_lon*dem_width;
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
    
    cal0=mean(lonlat)'; 
    xy_ratio=llh2local([cal0;0],[cal0+[0.1;0.1];0]);
    aspect_ratio=xy_ratio(1)/xy_ratio(2);

    y_posting=plot_pixel_m/plot_pixel_size*9e-6; %in degrees
    x_posting=y_posting/aspect_ratio;
    %plot_pixel_size=1;
    %plot_pixel_size=round(plot_pixel_m/plot_posting);

    
    
    
    % include the lon range in case user wants a larger area
    if ~isempty(lon_rg)
        x=[min(lon_rg)-2*x_posting:x_posting:max(lon_rg)+2*x_posting];
    else
        x=[min(lonlat(:,1))-2*x_posting:x_posting:max(lonlat(:,1))+2*x_posting];
    end
    % include the lon range in case user wants a larger area
    if ~isempty(lat_rg)
        y=[min(lat_rg)-2*y_posting:y_posting:max(lat_rg)+2*y_posting];
    else
        y=[min(lonlat(:,2))-2*y_posting:y_posting:max(lonlat(:,2))+2*y_posting];
    end
    
    [X,Y]=meshgrid(x,y);
    if bg_flag==0
        c=[[0 0 0];c]; % black background
    else
        c=[[1 1 1];c]; % white background
    end
    
    demx=round((lonlat(:,1)-x(1))/x_posting)+1;
    demy=round((lonlat(:,2)-y(1))/y_posting)+1;
    
    [~,uix]=unique([demy,demx],'rows');
    
    demx=demx(uix);
    demy=demy(uix);
    in_ph=in_ph(uix);
    col_ix=col_ix(uix);
    
    
    pixel_margin1=floor((plot_pixel_size-1)/2);
    pixel_margin2=ceil((plot_pixel_size-1)/2);
    
    ix1b=demy-pixel_margin1;
    ix1b(ix1b<1)=1;
    ix1e=demy+pixel_margin2;
    ix1e(ix1e>size(X,1))=size(X,1);     % DB
    ix2b=demx-pixel_margin1;
    ix2b(ix2b<1)=1;
    ix2e=demx+pixel_margin2;
    ix2e(ix2e>size(X,2))=size(X,2);     % DB

    R=zeros(size(X));
    
    
    nnix=~isnan(col_ix);
    col_ix= col_ix(nnix);
    in_ph=in_ph(nnix);
    demy=demy(nnix);       
    demx=demx(nnix);
    ix1b=ix1b(nnix); % DB    Bug fix 01/2017
    ix1e=ix1e(nnix); % DB    Bug fix 01/2017
    ix2e=ix2e(nnix); % DB    Bug fix 01/2017
    ix2b=ix2b(nnix); % DB    Bug fix 01/2017
    if plot_pixel_size==1
   	  for i=1 : length(in_ph)
          R(demy(i),demx(i))=col_ix(i)+1;
      end
    else
        for i=1 : length(in_ph)
    %        if ~(isnan(col_ix(i)))
    %             ix1=demy(i)-pixel_margin1:demy(i)+pixel_margin2;
    %             ix2=demx(i)-pixel_margin1:demx(i)+pixel_margin2;
    %             ix1=ix1(ix1>0&ix1<=size(X,1));
    %             ix2=ix2(ix2>0&ix2<=size(X,2));
    %             R(ix1,ix2)=col_ix(i)+1;


                  R(ix1b(i):ix1e(i),ix2b(i):ix2e(i))=col_ix(i)+1;
    %        end
        end
    end
    
    
    % optional to do extra filtering
    if strcmpi(filter_extra,'y')
        if strcmpi(filter_type,'median')

            R(R==0)=NaN;
            R2 = nanmedfilt2(R, [3 3]);
            R2(isnan(R)) = 0;
            R = R2;
            
            
        elseif strcmpi(filter_type,'mean')
            f = @(A)mean(A(~isnan(A)));
            R2 = nlfilter(R, [3 3], f);
            R2(isnan(R)) = 0;
            R = R2;
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
    if ref_radius<inf & ref_radius~=-inf
        p=plot(ref_centre(1),ref_centre(2),'k*');
        set(p,'markersize',10,'linewidth',1)
    end
    
   
    

elseif bg_flag==6 | bg_flag==7     % xy axes
    
    x_posting=plot_pixel_m;
    y_posting=plot_pixel_m;
    plot_pixel_size=1;

    x=[min(ps.xy(:,2))-2*x_posting:x_posting:max(ps.xy(:,2))+2*x_posting];
    y=[min(ps.xy(:,3))-2*y_posting:y_posting:max(ps.xy(:,3))+2*y_posting];
    
    [X,Y]=meshgrid(x,y);
    if bg_flag==7
        c=[[0 0 0];c]; % black background
    else
        c=[[1 1 1];c]; % white background
    end
    
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
    
 	for i=1:length(in_ph)
        if ~(isnan(col_ix(i)))
            p=plot(ps.xy(i,2),ps.xy(i,3),'.');
		    set(p,'color',c(col_ix(i),:));   
        end
    end
    axis equal
    axis tight
end
ph_lims=[max_ph,min_ph];







%%%% this is the part where I plot the external data.



% plotting of external data when requested
if plot_ext_data==1
    
    % Checking if there is a ph_disp variable stored in the files 
    if isfield(ext_data,{'ph_disp'})
        
        
        % plotting phase data
        col_ix2=round(((ext_data.ph_disp(:,1)-min_ph)*63/ph_range)+1);

        % accounting for values that are outside the extremes of the insar
        ix1 = find(col_ix2<2);
        ix2 = find(col_ix2>size(c_original,1));
        col_ix2(ix1) = 2;
        col_ix2(ix2) = size(c_original,1);
        ix = [ix1 ; ix2];
        ix_good = [1:size(col_ix2,1)]';
        ix_good(ix)=[];

        % plotting those values that are saturated by o marker
        if ~isempty(ix)
            p=scatter3(ext_data.lonlat(ix,1),ext_data.lonlat(ix,2),ext_data.ph_disp(ix,1),40,c_original(col_ix2(ix),:),'filled','o');
            set(p,'MarkerEdgeColor','w')
        end
        % plotting non-saturated values by square markers
        p=scatter3(ext_data.lonlat(ix_good,1),ext_data.lonlat(ix_good,2),ext_data.ph_disp(ix_good,1),40,c_original(col_ix2(ix_good),:),'filled','sq');
        set(p,'MarkerEdgeColor','w')
    else
        % plotting station data
        p=scatter(ext_data.lonlat(ix_good,1),ext_data.lonlat(ix_good,2),13,'^');
        set(p,'MarkerEdgeColor','w','MarkerfaceColor','k') 
        
    end
end

hold off
colormap(c);


