function []=coreg_pos(firstL,lastL,firstP,lastP,transL,transP,slc_osf)
% COREG_POS create a position file for fine coregistration
%    []=coreg_pos(firstL,lastL,firstP,lastP,transL,transP)
%
%   Andy Hooper, Aug 2005
%
% ======================================================================
% 07/2006 AH:  changed spacing of windows to give same number (3600) 
%              independent of image size
% 11/2009 AH:  code to use translations put back in 
% 03/2010 AH:  revert to just 3600 windows
% 09/2009 MA:  minor update for oversampled data + translation_lp
% 10/2010 JCM: drop negative values at the edges of the scene
% ======================================================================
%

if nargin<7
 slc_osf=1 ;              % default oversampling factor is 1
end
slc_osf
winSize=64*slc_osf ;      % correlation window size default: 64 should match with coreg.dorisin
                          %                                : 128 for oversampled data

nWinx=30    % # of windows in range
nWiny=120   % # of windows in azimuth


x=round(linspace(firstP+winSize+transP,lastP-winSize+transP,nWinx));
y=round(linspace(firstL+winSize+transL,lastL-winSize+transL,nWiny));

% drop negative values at the edges of the scene JCM
keep_ix=x>0;
keep_iy=y>0;
x=x(keep_ix);
y=y(keep_iy);


[Y,X]=meshgrid(y,x);

fid=fopen('fc_pos.in','w');
fprintf(fid,'%i  %i\n',[Y(:),X(:)]');
fclose(fid);

