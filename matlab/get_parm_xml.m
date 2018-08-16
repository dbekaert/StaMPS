function [parm_key,parm_key_comp] = get_parm_xml(xmlfile,parm)
% function which will read xml file and find the parameter_field matching
% the requested parameter. Will return string and numbers depending on the
% parameter_field
%
% Bekaert David - Jet Propulsion Laboratory
% modifications:
% DB    8/10/2016   Make sure not to end into an inf loop for the second while statement
% DB    3/11/2017   Trying to cpature incase there are mutiple variables with same parm matching

if nargin<2
    error('Require two inputs')
end

% initialize the field
counter=1;

% loop over the xml file 
fid = fopen(xmlfile,'r');
temp=1;
parm_key = 'FAIL';
parm_key_list = [];
parm_key_comp_list = [];
parm_key_comp = [];


while temp~=-1
    temp = fgetl(fid);
    if temp==-1
        break
    end
    
    
    % track the component to see multiple hits
    ix_comp =  strfind(temp,['<component']);
    if ~isempty(ix_comp)
        ix_comp =  strfind(temp,['=']);
        component = temp(ix_comp+1:end-1);
    end
    ix_comp_end =  strfind(temp,['</component']);
    if ~isempty(ix_comp_end)
        component = [];
    end
     
    % search for the variable
    ix =  strfind(temp,['"' parm '"']);
    if isempty(ix)
        ix = strfind(temp,['"' lower(parm) '"']);
    end
    if ~isempty(ix)
        % see if the value is at this string or not
        ix_value = strfind(temp,'<value>');
        while isempty(ix_value);
            temp = fgetl(fid);
            ix_value = strfind(temp,'<value>');
            if temp==-1
                break
            end
        end
        
        if temp~=-1
            % remove the start of <value>
            temp(1:ix_value+length('<value>')-1)=[];
            % remove the end of value
            ix_value = strfind(temp,'</value>');
            temp(ix_value:end)=[];
            parm_key_list{counter} = strtrim(temp);
            parm_key_comp_list{counter} = component;
            counter = counter+1;
        end
    end
end
fclose(fid);


% convert the array to a vector
if ~isempty(parm_key_list)
    clear parm_key
    for k=1:length(parm_key_list)
        if ~isempty(str2num(parm_key_list{k}))
            parm_key(k,1) = str2num(parm_key_list{k});
        else
            parm_key(k,:) = parm_key_list{k};
        end
    end
    try
        parm_key_comp = char(parm_key_comp_list);
    catch
        parm_key_comp = parm_key_comp_list;
    end
    
end

% see if the search was succesful
if strcmpi(parm_key,'FAIL')
    fprintf(['Could not find ' parm '... \n'])
end