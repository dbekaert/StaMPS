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
% ======================================================================


x=round(linspace(firstP+64,lastP-64,30));
y=round(linspace(firstL+64,lastL-64,120));
x=[firstP+64+transP:32:lastP-64+transP];
y=[firstL+64+transL:32:lastL-64+transL];
[Y,X]=meshgrid(y,x);

fid=fopen('fc_pos.in','w');
fprintf(fid,'%i  %i\n',[Y(:),X(:)]');
fclose(fid);
