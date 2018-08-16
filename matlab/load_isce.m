function [data_out] = load_isce(datafile)
% function which reads ISCE data, or data with .xml file into matlab.
% multiple bands are put in the thrid dimension.
% By Bekaert David - Jet Propulsion Laboratory
% modifications:
% 08/10/2016     DB  Change from vargout to 3D cube.
% 14/01/2017     DB  Include ISCE ifg support.
% 01/02/2017     DB  Include BIP for bands
% 23/04/2017     DB  Include Byte support to read mask files
% 04/11/2017     DB  Adding short support for DEM files

% the xml file
datafile_xml = [datafile '.xml'];
% check if the xml file exists
if exist(datafile_xml,'file')~=2
    error([datafile_xml ' does not exist'])
end

% getting relevant information 
[width] = get_parm_xml(datafile_xml,'width');
[length] = get_parm_xml(datafile_xml,'length');
[scheme] = get_parm_xml(datafile_xml,'scheme');
[number_bands] = get_parm_xml(datafile_xml,'number_bands');
[type] = get_parm_xml(datafile_xml,'data_type');


%% Different data cases
% check if this is an interferogram
if strcmpi(type,'CFLOAT')
    type_str = 'float32';
elseif strcmpi(type,'FLOAT')
    type_str = 'float32';
elseif strcmpi(type,'double')
    type_str = 'double';
elseif strcmpi(type,'byte')
    type_str = 'uint8';
elseif strcmpi(type,'short')
    type_str = 'short';
else
    error([type ' to be included...'])
end

if strcmpi(type,'CFLOAT')
    if strcmpi(scheme,'BSQ')
        if number_bands==1
            fid = fopen(datafile,'r');
            data_out = fread(fid,[2*width inf],type_str);
            fclose(fid);
            data_out = complex(data_out(1:2:end,:),data_out(2:2:end,:));
        else
            error('To be coded further...')
        end
    elseif strcmpi(scheme,'BIP')
        if number_bands==1
            fid = fopen(datafile,'r');
            data_out = fread(fid,[2*width inf],type_str);
            fclose(fid);
            data_out = complex(data_out(1:2:end,:),data_out(2:2:end,:));
        else
            error('To be coded further...')
        end
    else          
        error([type ': ' scheme ' to be included...'])
    end   
else
    
    if strcmpi(scheme,'BIP')
        if number_bands==1
            fid = fopen(datafile,'r');
            data_out = fread(fid,[width inf],type_str);
            
            clear data
            fclose(fid);
        else
            fid = fopen(datafile,'r');
            data = fread(fid,[number_bands*width inf],type_str);
            fclose(fid);
            for k=1:number_bands
                data_out(:,:,k) = data(k:number_bands:end,:);
            end
        end
    elseif strcmpi(scheme,'BIL')
        fid = fopen(datafile,'r');
        data = fread(fid,[width inf],type_str);
        fclose(fid);
        for k=1:number_bands
            data_out(:,:,k) = data(:,k:number_bands:end);
        end
        clear data
    else
        error([type ': ' scheme ' to be included...'])
    end
end

