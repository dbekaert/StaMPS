function [values]=readparm(fname,parm,numval,log_flag)
%READPARM read a parameter from a file
%   [values]=readparm(fname,parm,numval,log_flag)
%   Andy Hooper, June 2011
%
%   ======================================================================
%   09/2006 AH: 'v' option added

%
if nargin<3
    numval=1;
end

if nargin<4
    log_flag=1;
end

fid=fopen(fname);
parms=textscan(fid,'%s');
fclose(fid);

parms=parms{1};
ix=strmatch(parm,parms);
msg=[parms{ix},'='];
for i=1:numval
    parmval=parms{ix+i};
    msg=[msg,parmval,' '];
    if isempty(str2num(parmval))
        values{i}=parmval;
    else        
        values{i}=str2num(parmval);
    end
end

if log_flag
    logit(msg)
end

if numval==1
    values=values{1};
end
