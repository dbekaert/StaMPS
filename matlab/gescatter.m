function gescatter(filename,lon,lat,c,varargin)
% GESCATTER - create a scatter plot in Google Earth
%
%   GESCATTER(FILENAME,LON,LAT,C) - creates a .kml file that
%   displays colored circles at the locations specified by the 
%   vectors LON and LAT similar to ML's builtin function, SCATTER.
%   The color of the circles is scaled relative to the
%   values provided in third input, C.
%
%   OPTIONS AND SYNTAX - Optional inputs are entered as 
%   property/value pairs.  Valid optional properties are:
%
%   GESCATTER(...,'colormap','hot') - uses Matlabs 'hot'
%   colormap instead of the default (jet). Also accepts function
%   handles (@hot), or custom colormaps (m x 3 matrices).
%
%   GESCATTER(...,'clims',[low high]) - limit the color
%   values to a specified values (similar to CAXIS on a ML
%   figure). Clims should be supplied as a 2-element array.
%
%   GESCATTER(...,'time',timevector) - assigns a time to each
%   point. The length of the timevector array should be the same
%   as LAT, LON, and C.
%
%   GESCATTER(...,'scale',size) - scales the size of the dots in the
%   Google Earth file.  Default value is 0.4.
%
%   EXAMPLE
%
%   %generate some data
%   x=(0:0.05:6*pi);
%   lon = -122.170087 + cos(x)*0.01;
%   lat = 37.455697 + x*0.001;
%
%   %color the points according to their latitude
%   gescatter('foo.kml',lon,lat,lat)
%
% SEE ALSO scatter

% A. Stevens @ USGS 3/04/2009
% astevens@usgs.gov
% 
% Modifications:
% M. Arikan @ TU Delft  - OPACITY option added 10/04/2010
% A Hooper - minor changes 17 March 2015 

% Copyright (c) 2010, Andrew Stevens
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.


%default values
clims=[min(c) max(c)];
cmap=fliplr(jet);
t=[];
scale=0.2;
opacity=0.4;  % default value for 40 percent opacity [MA]

%parse inputs and do some error-checking
if nargin>0
    [m,n]=size(varargin);
    opts={'clims','time','scale','colormap','opacity'};
    for i=1:n;
        indi=strcmpi(varargin{i},opts);
        ind=find(indi==1);
        if isempty(ind)~=1
            switch ind
                case 1
                    clims=varargin{i+1};
                    if numel(clims)~=2
                        error('Clims should be a two-element array.')
                    end
                    
                case 2
                    t=varargin{i+1};
                    if any(isnan(t))
                        error('Time vector should not contain NaNs.')
                    end
                    if ~isnumeric(t)
                        error('Time should be entered in ML Datenum format.')
                    end
                case 3
                    scale=varargin{i+1};
                case 4
                    cmap=varargin{i+1};
                    %check size of numeric colormap input
                    if isa(cmap,'numeric')
                        [m,n]=size(cmap);
                        if n~=3
                            error('Custom colormap must have 3 columns.')
                        end
                        cmap=fliplr(cmap);
                    else
                        %if standard colormap is supplied
                        if isa(cmap,'function_handle')
                            cmap= func2str(cmap); 
                        end
                        cmap=fliplr(feval(cmap));
                    end
                case 5
                    opacity=varargin{i+1};
                    if ~isnumeric(opacity) ||  opacity > 1.0 || opacity < 0
                        error('Opacity should be numeric in 0.0 - 1.0 range.')
                    end
                    
            end
        end
    end
end


[pathstr,namer] = fileparts(filename);

%get rid on nans
gind=(isfinite(lon) & isfinite(lat) & isfinite(c));
lon=lon(gind);
lat=lat(gind);
c=c(gind);


%figure out the rgb colors of each value
cvals=[-inf;linspace(clims(1),clims(2),...
    length(cmap(:,1))-2)';inf];
[n,bin]=histc(c,cvals);
colors=cmap(bin,:);

% OPACITY
opacity=repmat(opacity,size(colors,1),1);
colors=[opacity colors];

%convert to GE's hex format
rgb=cellfun(@(x)(dec2hex(floor(x.*255),2)),...
    num2cell(colors),'uni',0);

%write the GE file
header=['<?xml version="1.0" encoding="UTF-8"?>',...
    '<kml xmlns="http://www.opengis.net/kml/2.2">',...
    '<Document><name>',namer,'</name>'];
footer='</Document></kml>';

h = waitbar(0,'Creating file, Please wait...');
set(h,'name','Creating Google Earth file')

fid = fopen(filename, 'wt');
fprintf(fid, '%s \n',header);

for i=1:length(lon)
    
    %create a style to hide each point in one document
    fprintf(fid,'%s \n','<Style id="folderStyle">');
    fprintf(fid,'%s \n','<ListStyle>');
    fprintf(fid,'%s \n','<listItemType>checkHideChildren</listItemType>');
    fprintf(fid,'%s \n','</ListStyle>');
    fprintf(fid,'%s \n','</Style>');
    
    %define the point style
    fprintf(fid,'%s \n','<Style id="cpoint">');
    fprintf(fid,'%s \n','<IconStyle>');
    %fprintf(fid,'%s \n',['<color>ff',[rgb{i,:}],'</color>']);
    fprintf(fid,'%s \n',['<color>',[rgb{i,:}],'</color>']);
    fprintf(fid,'%s \n',['<scale>',sprintf('%.1f',scale),'</scale>']);
    fprintf(fid,'%s \n',['<Icon><href>http://maps.google.com/mapfiles/',...
        'kml/shapes/shaded_dot.png</href></Icon>']);
    fprintf(fid,'%s \n','</IconStyle>');
    fprintf(fid,'%s \n','</Style>');
    
    %add the placemark
    fprintf(fid, '%s \n','<Placemark>');
    fprintf(fid,'%s \n','<styleUrl>#cpoint</styleUrl>');
    
    %create a simple description for each point
    fprintf(fid, '%s \n','<description><![CDATA[<table width="200"></table>');
    fprintf(fid, '%s \n',['<h2>Filename: ',namer,'<br>']);
    fprintf(fid, '%s \n',['<h3>Value: ',sprintf('%.1f',c(i)),'<br>']);
    if ~isempty(t)
        fprintf(fid, '%s \n',['Time (GMT): ',datestr(t(i)),'<br>']);
    end
    fprintf(fid, '%s \n',']]></description>');
    
    
    fprintf(fid,'%s \n','<Point>');
    fprintf(fid,'%s','<coordinates>');
    fprintf(fid, ' %.6f, %.6f, %.2f', [lon(i) lat(i) c(i)]);
    fprintf(fid,'%s \n','</coordinates>');
    fprintf(fid,'%s \n','</Point>');
    
    if ~isempty(t)
        fprintf(fid,'%s \n','<TimeSpan>');
        fprintf(fid,'%s \n',['<begin>',datestr(t(1),29),...
            'T',datestr(t(1),13),'Z</begin>']);
        fprintf(fid,'%s \n',['<end>',datestr(t(end),29),...
            'T',datestr(t(end),13),'Z</end>']);
        fprintf(fid,'%s \n','</TimeSpan>');
    end
    
    
    fprintf(fid, '%s \n','</Placemark>');
    
    %waitbar(i/length(lon),h,sprintf('%d%% complete...',...
        %round((i/length(lon))*100)));
    
end

fprintf(fid, '%s \n','<styleUrl>#folderStyle</styleUrl>');
fprintf(fid, '%s \n',footer);

close(h);
fclose(fid);
