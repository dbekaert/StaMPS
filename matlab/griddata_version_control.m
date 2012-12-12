function [z] = griddata_version_control(gridX,gridY,gridZ,x,y,method,matlab_version)
% This function just computes the gridded data set, but its input is
% modified to work with matlab version newer than 2012.
% By David Bekaert - December 2012


if matlab_version>=2012
    z=griddata(gridX,gridY,gridZ,x,y,method);
else
    z=griddata(gridX,gridY,gridZ,x,y,method,{'QJ'});
end
