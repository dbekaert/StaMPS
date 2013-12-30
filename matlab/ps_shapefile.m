function ps_shapefile(filename,data,step)
%
%   12/2013 MA: Initial version
%
%

if nargin < 3
  step=':';
end

% get velocity
%if ~exist('ps_plot_v-d.mat')
% ps_plot('v-d',-1)
%end
%load ps_plot_v-d  ph_disp 
%data= ph_disp;

disp('Loading lonlat matrix...')
load ps2 lonlat


step
if isinteger(step)
 step=1:step:size(lonlat,1);
end

lon=lonlat(step,:,1);
lat=lonlat(step,:,2);
data=double(data(step,:));

shape=struct('Geometry', 'Point','X',num2cell(lon),'Y',num2cell(lat),'V',num2cell(V));

shapewrite(shape,filename);