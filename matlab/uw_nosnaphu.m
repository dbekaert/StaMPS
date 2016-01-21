function [ph_uw]=uw_nospnaphu(ph,xy,day,options)
%UW_NOSNAPHU unwrap first in time then space
%   UW_NOSNAPHU(PH,XY,DAY,ALPHA,MAX_RM_FRACTION)
%
%   Based on UW_PSEUDO3D by Andy Hooper from June 2006
%
%   Andy Hooper, Sep 2015

uw_triangulate(ph,xy,day)
uw_unwrap_time(options.time_win)
ph_uw=uw_unwrap_space;