function [ref_ps]=ps_setref()
%PS_SETREF find reference PS
%
%   Andy Hooper, June 2006
%
%   07/2006 AH changed to use reference lon/lat if set
%   02/2010 AH option to also use circular reference area


load psver
psname=['ps',num2str(psver)];

ps2=load(psname);

[ref_lon,parmname]=getparm('ref_x');
if strcmp(parmname,'ref_x')
    ref_x=getparm('ref_x');
    ref_y=getparm('ref_y');
    ref_ps=find(ps2.xy(:,2)>ref_x(1)&ps2.xy(:,2)<ref_x(2)&ps2.xy(:,3)>ref_y(1)&ps2.xy(:,3)<ref_y(2));
else
    ref_lon=getparm('ref_lon');
    ref_lat=getparm('ref_lat');
    ref_centre_lonlat=getparm('ref_centre_lonlat');
    ref_radius=getparm('ref_radius');
    ref_ps=find(ps2.lonlat(:,1)>ref_lon(1)&ps2.lonlat(:,1)<ref_lon(2)&ps2.lonlat(:,2)>ref_lat(1)&ps2.lonlat(:,2)<ref_lat(2));
    if ref_radius<inf
        ref_xy=llh2local(ref_centre_lonlat',ps2.ll0)*1000;
        xy=llh2local(ps2.lonlat(ref_ps,:)',ps2.ll0)*1000;
        dist_sq=(xy(1,:)-ref_xy(1)).^2+(xy(2,:)-ref_xy(2)).^2; 
        ref_ps=ref_ps(dist_sq<=ref_radius^2);
    end
end

if isempty(ref_ps)
   ref_ps=[1:ps2.n_ps]';
end

disp([num2str(length(ref_ps)),' ref PS selected'])





