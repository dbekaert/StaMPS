function ps_gescatter(filename,data,step,opacity)
% Generate kml from PS results
%
%
%
%   Andy Hooper, June 2010
%
%   ======================================================================
%   10/2010 MA: Initial version
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

step   % disseminate data
if nargin < 3
  step=1;
end
step   % disseminate data

disp('Loading lotlat matrix...')
load ps2 lonlat

step=1:step:size(lonlat,1);

gescatter(filename,lonlat(step,1),lonlat(step,2),data(step,:),'size',5,'colormap','jet','opacity',opacity) 
 
%EOF
