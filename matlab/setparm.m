function []=setparm(parmname,value,newflag)
%SETPARM sets parameters in a saved workspace
%   SETPARM(PARMNAME,VALUE,FLAG) 
%   Only enough characters of PARMNAME to make it unique need be typed.
%   If VALUE is set to nan, the parameter is reset to the default value.   
%   FLAG is optional, valid values are:
%       1 = add a new parameter (PARMNAME must be typed in full)
%       2 = add a new parameter to the local parameter workspace
%      -1 = delete parameter from workspace(VALUE ignored)
%      -2 = delete parameter from local workspace (VALUE ignored)
%
%   Andy Hooper, June 2006
%
%   07/2006 AH Add patch compatibilty
%   03/2007 AH Support for local parameter workspace added
%   10/2007 AH Parameters displayed in alphabetical order
%   12/2007 AH Add option to reset to a parameter to the default value
%   03/2008 AH Default processing amended


parmfile='parms';
localparmfile='localparms.mat';

if exist('./parms.mat','file')
    parms=load(parmfile);
elseif exist('../parms.mat','file')
    parmfile='../parms';
    parms=load(parmfile);
else
    parms=struct();
end

if exist(localparmfile,'file')
    localparms=load(localparmfile);
else
    localparms=struct('Created',date);
end

if nargin>2 & newflag==1
    if isnumeric(value)
        disp([parmname,' = ',num2str(value)])
    else
        disp([parmname,' = ',value])
    end
    parms=setfield(parms,parmname,value);
    save(parmfile,'-struct','parms')
      
else
    if nargin>1
        if ~isnumeric(parmname)
            parmnum=strmatch(parmname,fieldnames(parms)); 
            if length(parmnum)>1
                error(['Parameter ',parmname,'* is not unique'])
            elseif isempty(parmnum)
                error(['Parameter ',parmname,'* does not exist'])
            end
        else
            parmnum=parmname;
        end
        parmnames=fieldnames(parms);
        if size(parmnames,1)<parmnum
            error(['There are only ',num2str(size(parmnames,1)),' fields'])
        end
        parmname=parmnames{parmnum};
        if isnan(value)
            if strcmpi(parmname,'small_baseline_flag')
                error('Default reset option not possible for small_baselines_flag')
            end
            parms=rmfield(parms,parmname);
            save(parmfile,'-struct','parms')
            if isfield(localparms,parmname)
                localparms=rmfield(localparms,parmname);
                save(localparmfile,'-struct','localparms')
            end
            % AH 03/2008
            %if strcmpi(parms.small_baseline_flag,'y') 
            %    sb_parms_default
            %end
            ps_parms_default
            value=getparm(parmname);
            disp([parmname,' reset to default value'])
        end
        if isnumeric(value)
            disp([parmname,' = ',num2str(value)])
        else
            disp([parmname,' = ',value])
        end
        if nargin>2 
            if newflag==2
                localparms=setfield(localparms,parmname,value);
                save(localparmfile,'-struct','localparms')
                disp('Added to LOCAL parameter file')
            elseif newflag==-1
                parms=rmfield(parms,parmname);
                save(parmfile,'-struct','parms')
                disp([parmname,' removed from parameter file'])
            elseif newflag==-2
                localparms=rmfield(localparms,parmname);
                if size(fieldnames(localparms),1)>1
                    save(localparmfile,'-struct','localparms')
                else
                    delete(localparmfile)
                end
                disp([parmname,' removed from LOCAL parameter file'])
                
            else
                error('Invalid value for NEWFLAG')
            end
        else
            if ~isfield(localparms,parmname)
                parms=setfield(parms,parmname,value);
                save(parmfile,'-struct','parms')
            else
                localparms=setfield(localparms,parmname,value);
                save(localparmfile,'-struct','localparms')
                disp('Warning: Only LOCAL parameter file updated')
            end
        end
    elseif nargin >0
        error('Format is: SETPARM(PARMNAME,VALUE,[NEWFLAG])')
    else
        disp(orderfields(parms))
        if size(fieldnames(localparms),1)>1
            localparms
        end
    end
end


    
    

