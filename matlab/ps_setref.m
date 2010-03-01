function [ref_ps]=ps_setref()
%PS_SETREF find reference PS
%
%   Andy Hooper, June 2006
%
%   07/2006 AH changed to use reference lon/lat if set


load psver
psname=['ps',num2str(psver)];

ps2=load(psname);

[ref_lon,parmname]=getparm('ref_lon');
if ~strcmp(parmname,'ref_lon')
    ref_x=getparm('ref_x');
    ref_y=getparm('ref_y');
    ref_ps=find(ps2.xy(:,2)>ref_x(1)&ps2.xy(:,2)<ref_x(2)&ps2.xy(:,3)>ref_y(1)&ps2.xy(:,3)<ref_y(2));
else
    ref_lat=getparm('ref_lat');
    ref_ps=find(ps2.lonlat(:,1)>ref_lon(1)&ps2.lonlat(:,1)<ref_lon(2)&ps2.lonlat(:,2)>ref_lat(1)&ps2.lonlat(:,2)<ref_lat(2));
end

if isempty(ref_ps)
   ref_ps=[1:ps2.n_ps]';
end

disp([num2str(length(ref_ps)),' ref PS selected'])





