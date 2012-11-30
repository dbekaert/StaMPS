function []=logit(logmsg,whereto,parent_flag)
%LOGIT write message to log file and/or stdout
%   [RC]=LOGIT(LOGMSG,LEVEL)
%   LEVEL = 0 (default) writes to stamps.log and stdout
%   LEVEL = 1 writes to stamps.log only
%   LEVEL = 2 writes to stdout only
%   LEVEL = 3 writes to debug.log only
%   LEVEL = FILENAME writes to FILENAME and stdout
%   PARENT_FLAG = 0 (default) writes to current directory
%   PARENT_FLAG = 1 writes to parent directory
%
% Andy Hooper, March 2010
%
%   ================================================================
%   11/2012 AH: Add option to specify log filename
%   ================================================================
   
[fstack]=dbstack ;
if size(fstack,1)>1
    fname=upper(fstack(2).name);
else
    fname='Command line';
end

if nargin<1
    logmsg=0;
end

if nargin<2
    whereto=0;
end

if nargin<3
    parent_flag=0;
end

if isnumeric(logmsg)
   switch logmsg
   case 0
       logmsg='Starting';
   case 1
       logmsg='Finished';
   end
end

if strcmp(logmsg(end-1:end),'\n')
    logmsg=logmsg(1:end-2);
end

if ~isnumeric(whereto)
    logfile=whereto;
    whereto=0;
else
    logfile='STAMPS.log';
end

debugfile='DEBUG.log';
if parent_flag==1
   logfile=['../',logfile];
   debugfile=['../',debugfile];
end

if whereto<2
    logid=fopen(logfile,'a');
    if logid>0
        fprintf(logid,'%s %-16s %s\n',datestr(now),fname,logmsg);
        fclose(logid);
    end
end

if whereto==0 | whereto==2
    fprintf('%s: %s\n',fname,logmsg);
end

if whereto==3
    logid=fopen(debugfile,'a');
    if logid>0
        fprintf(logid,'%s %-16s %s\n',datestr(now),fname,logmsg);
        fclose(logid);
    end
end

