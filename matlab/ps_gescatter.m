function ps_gescatter(filename,data,step)
% Generate kml form PS results
%
%
%
%   Andy Hooper, June 2006
%
%   ======================================================================
%   10/2010 MA: Initial version
%   ======================================================================

% TODO
% 1. quadtree zoom levels
% 2. attach TS series for pixels
% 3. add colorbar

% get velocity
%if ~exist('ps_plot_v-d.mat')
% ps_plot('v-d',-1)
%end
%load ps_plot_v-d  ph_disp 
%data= ph_disp;

if nargin < 3
  step=1;
end
step   % disseminate data

load ps2 lonlat

step=1:step:size(lonlat,1);

gescatter(filename,lonlat(step,1),lonlat(step,2),data(step,:),'size',5,'colormap','jet','opacity',0.4) 
 
%EOF
