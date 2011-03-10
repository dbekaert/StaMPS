function putdates(xstart, ystart, labels, labeloffset, fontsize)
% putdate to ts_plot for PS_PLOT function
%
%
%
%   Andy Hooper, June 2006
%
%   ======================================================================
%   11/2010 MA: initial version for StAMPS time series 
%   03/2011 AH: fix bug
%   ======================================================================
     axis off
     rectangle('Facecolor','w')
     %labels=dates(ind_Btemp);
     posx = xstart; % 0.05;
     posy = ystart; %1;
     %labeloffset= 0.035;  % 0.025 decrease this to fit more labels % label line offset
                           % text font 8   labelloffset: 0.03 or 0.025
                           %           9                 0.035
     for i=1:size(labels,1)
         posy= posy - labeloffset;
         if posy < 0.01
             posy = 1 - labeloffset;  % initialize again for next label coln
             posx = posx + 0.5 ;
         end
         %text(posx, posy, [  num2str(i,'%2d') '. ' char(labels(i)) ],...%c
         text(posx, posy, [  num2str(i,'%2d') '. ' char(labels(i,:)) ],...
             'Units','normalized','Fontsize',fontsize)
     end
