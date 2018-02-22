function ps_gescatter(filename,data,step,opacity,climsmtx,lon_rg,lat_rg,stdev_cutoff)
% Generate kml from PS results
%
%
%
%   Andy Hooper, June 2010
%
%   ======================================================================
%   10/2010 MA: Initial version
%   01/2011 MA: define color limits
%   06/2017 DB: Take random selection instead of fixed step-size to avoid patterns
%   ======================================================================

% TODO list
% 1. quadtree zoom levels
% 2. attach TS series for pixels
% 3. add colorbar

% get velocity
%if ~exist('ps_plot_v-d.mat')
% ps_plot('v-d',-1)
%end
%load ps_plot_v-d  ph_disp 
%data= ph_disp;

if nargin < 4
  opacity=0.4;
end

if nargin < 3
  step=1;
end
if nargin<5
    climsmtx=[];
end
if nargin<6 | isempty(lon_rg)
    lon_rg=[-inf,inf];
end
if nargin<7 | isempty(lat_rg)
    lat_rg=[-inf,inf];
end
if nargin<8 | isempty(stdev_cutoff)
    stdev_cutoff=inf;
end
step   % print increments of disseminated data

disp('Loading lonlat matrix...')
load ps2 lonlat n_ps

  


ix=lonlat(:,1)>=min(lon_rg)&lonlat(:,1)<=max(lon_rg)&lonlat(:,2)>=min(lat_rg)&lonlat(:,2)<=max(lat_rg);
if stdev_cutoff<inf
    load ./mv2.mat mean_v_std
    ix2=mean_v_std<stdev_cutoff;
    ix=ix&ix2;
end
lonlat=lonlat(ix,:);
if size(data,1)==1
    data=data.';
end
data=data(ix,:);

A = randperm(n_ps);
step = sort(A(1:ceil(n_ps./step)));
clear A


%step=1:step:size(lonlat,1);

%gescatter(filename,lonlat(step,1),lonlat(step,2),data(step,:),'size',5,'colormap','jet','opacity',opacity) 

%jetme=fliplr(jet);  % change colors order for climsmtx
aa=jet(63);
jetme=flipud([repmat(aa(8,:),15,1);repmat(aa(23,:),10,1);repmat([0.3,1,0.3],14,1);repmat(aa(40,:),10,1);repmat(aa(56,:),15,1)]);

% size value range for google icon (shaded_dot): 0.0 to 1.0

if ~isempty(climsmtx) % do clims
  gescatter(filename,lonlat(step,1),lonlat(step,2),data(step,:),'scale',0.4,'colormap',jetme,'opacity',opacity,'clims',climsmtx) 
else            % don't do clims
  gescatter(filename,lonlat(step,1),lonlat(step,2),data(step,:),'scale',0.2,'colormap','jet','opacity',opacity)
end
 
%EOF
