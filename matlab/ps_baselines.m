function [bperp,ix]=ps_baselines(max_bperp)
%PS_BASELINES returns perpendicular baselines above a threshold
%  [BPERP,IX]=PS_BASELINES(MAX_BPERP)
%
%   Andy Hooper, July 2009
%
%   ======================================================================
%   11/2010 AH: Return baselines above instead of below threshold
%   ======================================================================


if nargin<1
    max_bperp=inf;
end

load psver
psname=['ps',num2str(psver)];
load(psname,'bperp');
ix=find(abs(bperp)>=max_bperp)';
bperp=bperp(ix);
