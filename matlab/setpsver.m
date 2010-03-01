function []=setpsver(new_psver)
% SETPSVER sets psver
%
%   Andy Hooper, June 2006


load psver
disp(['psver currently: ',num2str(psver)])

if nargin>0
    psver=new_psver;
    save psver psver
    disp(['psver now set to: ',num2str(psver)])
end

