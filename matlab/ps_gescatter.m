function ps_gescatter(filename,data,step,opacity,climsmtx)
% Generate kml from PS results
%
%
%
%   Andy Hooper, June 2010
%
%   ======================================================================
%   10/2010 MA: Initial version
%   01/2011 MA: define color limits
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
step   % print increments of disseminated data

disp('Loading lonlat matrix...')
load ps2 lonlat

step=1:step:size(lonlat,1);

%gescatter(filename,lonlat(step,1),lonlat(step,2),data(step,:),'size',5,'colormap','jet','opacity',opacity) 

jetme=fliplr(jet);  % change colors order for climsmtx

if (nargin > 4) % do clims
  gescatter(filename,lonlat(step,1),lonlat(step,2),data(step,:),'size',5,'colormap',jetme,'opacity',opacity,'clims',climsmtx) 
else            % don't do clims
  gescatter(filename,lonlat(step,1),lonlat(step,2),data(step,:),'size',5,'colormap','jet','opacity',opacity)
end
 
%EOF
