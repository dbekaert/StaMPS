function [value,parmname]=getparm(parmname,printflag)
%GETPARM get parameter value from parms.mat
%   GETPARM(PARMNAME) 
%   Only enough characters of PARMNAME to make it unique need be typed
%
%   Andy Hooper, July 2006
%
%   10/2007 AH Parameters displayed in alphabetical order

ps_parms_default

if nargin<2
    printflag=0;
end

parmfile='parms';
localparmfile='localparms';

if exist('./parms.mat','file')
    parms=load(parmfile);
elseif exist('../parms.mat','file')
    parmfile='../parms';
    parms=load(parmfile);
else
    error('parms.mat not found')
end

if exist('localparms.mat','file')
    localparms=load(localparmfile);
else
    localparms=struct('Created',date);
end


if nargin < 1
    disp(orderfields(parms))
    if size(fieldnames(localparms),1)>1
        localparms
    end
else
    parmnum=strmatch(parmname,fieldnames(parms)); 
    if length(parmnum)>1
        error(['Parameter ',parmname,'* is not unique'])
    elseif isempty(parmnum)
        parmname=[];
        value=[];
    else
        parmnames=fieldnames(parms);
        parmname=parmnames{parmnum};
        if isfield(localparms,parmname)
            value=getfield(localparms,parmname);
        else
            value=getfield(parms,parmname);
        end
    end
    if printflag~=0
        if isnumeric(value)
            fprintf(['   PARM: %s=',repmat('%g ',1,200)],parmname,value)
            fprintf('\n')
        else
            fprintf('   PARM: %s=''%s''\n',parmname,value)
        end
    end
end



    
    

