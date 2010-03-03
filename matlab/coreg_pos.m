function []=coreg_pos(firstL,lastL,firstP,lastP,transL,transP)
% COREG_POS create a position file for fine coregistration
%    []=coreg_pos(firstL,lastL,firstP,lastP,transL,transP)
%
%   Andy Hooper, Aug 2005
%
% ======================================================================
% 07/2006 AH: changed spacing of windows to give same number (3600) 
%              independent of image size
% 11/2009 AH: code to use translations put back in 
% 03/2010 AH: revert to just 3600 windows
% ======================================================================

x=round(linspace(firstP+64+transP,lastP-64+transP,30));
y=round(linspace(firstL+64+transL,lastL-64+transL,120));
[Y,X]=meshgrid(y,x);

fid=fopen('fc_pos.in','w');
fprintf(fid,'%i  %i\n',[Y(:),X(:)]');
fclose(fid);

