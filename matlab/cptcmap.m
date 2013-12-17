function varargout = cptcmap(varargin)
%CPTCMAP Apply a .cpt file as colormap to an axis
%
% cptcmap(name);
% cptcmap(name, ax);
% cptcmap(... param, val, ...);
% [cmap, lims,ticks,bfncol] = cptcmap(...)
%
% This function creates and applies a colormap defined in a color palette
% table (.cpt file).  For a full description of the cpt file format, see
% the Generic Mapping Tools documentation (http://gmt.soest.hawaii.edu/).
% Color palette files provide more flexible colormapping than Matlab's
% default schemes, including both discrete and continuous gradients, as
% well as easier direct color mapping.
%
% Limitations: X11 color names not supported, patterns not supported, CMYK
% not supported yet
%
% Input variables:
%
%   name:       .cpt file name.  You may either specify either the full
%               path, or just the file name.  In the latter case, the
%               function will look for the file in the folder specified by
%               the cptpath variable in the first line of code; please
%               modify this variable to match your own cpt folder.
%
%   ax:         handle of axis or axes where colormap should be applied
%               (colormaps will effect the entire figure(s), but axis clim
%               adjustments for direct scaling will only affect the
%               specified axes).  If no axis is specified and no output
%               variables are supplied, colormap will be applied to the
%               current axis.  If no axis is specified and output variables
%               are supplied, the colormap will not be applied to any axes.
%
%   'showall':  When this option is used, a figure is created displaying
%               colorbars for all colormaps contained in the .cpt folder.
%
% Optional input variables (passed as parameter/value pairs):
%
%   'mapping':  'scaled' or 'direct'.  Scaled mapping spreads the colormap
%               to cover the color limits of the figure.  Direct mapping
%               resets the color limits of the axes so that colors are
%               mapped to the levels specified by the .cpt file. ['scaled']
%
%   'ncol':     number of colors in final colormap. If not included or NaN,
%               this function will try to choose the fewest number of
%               blocks needed to display the colormap as accurately as
%               possible. I have arbitrarily chosen that it will not try to
%               create more than 256 colors in the final colormap when
%               using this automatic scheme.  However, you can manually set
%               ncol higher if necessary to resolve all sharp breaks and
%               gradients in the colormap.
%
%   'flip':     if true, reverse the colormap order [false]
%
% Output variables:
%
%   cmap:       ncol x 3 colormap array
%
%   lims:       1 x 2 array holding minimum and maximum values for which
%               the colormap is defined.  
%
%   ticks:      vector of tick values specifying where colors were defined
%               in the original file
%
% Example:
%
%   [lat, lon, z] = satbath(10);
%   pcolor(lon, lat, z);
%   shading flat;
%   cptcmap('GMT_globe', 'mapping', 'direct');
%   colorbar; 

% Copyright 2010 Kelly Kearney
%
% Modifications:
% David Bekaert 	12/2013 	Put the cpt folder hardcoded at the same location as the script

%------------------------------
% Parse input
%------------------------------

% Set up .cpt file folder [DB]
[temp] = which('cptcmap.m');
[path, b, c] = fileparts(temp);
cptpath = [path filesep 'cptfiles' filesep];
if ~exist(cptpath, 'dir')
    error('Please modify the cptpath variable in this code to point to the directory where your .cpt files are stored');
end

% Check for 'showall' option first

if nargin == 1 & strcmp(varargin{1}, 'showall')
    plotcmaps(cptpath);
    return;
end

% Name of file

if exist(varargin{1}, 'file')   % full filename and path given
    filename = varargin{1};
else                            % only file name given
    if strcmp(varargin{1}(end-3:end), '.cpt')   % with extension
        filename = fullfile(cptpath, varargin{1});
    else                                        % without extension
        filename = fullfile(cptpath, [varargin{1} '.cpt']);   
    end
    if ~exist(filename, 'file')
        error('Specified .cpt file not found');
    end
end

% Axes to which colormap will be applied

if nargin > 1 && isnumeric(varargin{2}) && all(ishandle(varargin{2}(:)))
    ax = varargin{2};
    pv = varargin(3:end);
    applycmap = true;
elseif nargout == 0
    ax = gca;
    pv = varargin(2:end);
    applycmap = true;
else
    pv = varargin(2:end);
    applycmap = false;
end

% Optional paramter/value pairs
    
Opt = struct('mapping', 'scaled', ...
             'ncol', NaN, ...
             'flip', false);
         
Opt = parsepv(Opt, pv);
     
%------------------------------
% Calculate colormap and apply
%------------------------------

[cmap, lims,ticks,bfncol] = cpt2cmap(filename, Opt.ncol);
if Opt.flip
    if strcmp(Opt.mapping, 'direct')
        warning('Flipping colormap with direct mapping may lead to odd color breaks');
    end
    cmap = flipud(cmap);
end

if applycmap
    for iax = 1:numel(ax)
        axes(ax(iax));
        if strcmp(Opt.mapping, 'direct')
            set(ax(iax), 'clim', lims);
        end
        colormap(cmap);
    end
end

%------------------------------
% Output
%------------------------------

if nargout == 1
    varargout{1} = cmap;
elseif nargout == 2
    varargout{1} = cmap;
    varargout{2} = lims;
elseif nargout == 3
    varargout{1} = cmap;
    varargout{2} = lims;
    varargout{3} = ticks;
elseif nargout == 4
    varargout{1} = cmap;
    varargout{2} = lims;
    varargout{3} = ticks;
    varargout{4} = bfncol;
end

%------------------------------
% Subfunction: Read colormap 
% from file
%------------------------------

function [cmap, lims, ticks, bfncol] = cpt2cmap(file, ncol)

% Read file

fid = fopen(file);
txt = textscan(fid, '%s', 'delimiter', '\n');
txt = txt{1};
fclose(fid);

isheader = strncmp(txt, '#', 1);
isfooter = strncmp(txt, 'B', 1) | strncmp(txt, 'F', 1) | strncmp(txt, 'N', 1); 

% Extract color data, ignore labels (errors if other text found)

ctabletxt = txt(~isheader & ~isfooter);
ctable = str2num(strvcat(txt(~isheader & ~isfooter)));
if isempty(ctable)
    nr = size(ctabletxt,1);
    ctable = cell(nr,1);
    for ir = 1:nr
        ctable{ir} = str2num(strvcat(regexp(ctabletxt{ir}, '[\d\.-]*', 'match')))';
    end
    try 
        ctable = cell2mat(ctable);
    catch
        error('Cannot parse this format .cpt file yet');
    end 
end

% Determine which color model is used (RGB, HSV, CMYK, names, patterns,
% mixed)

[nr, nc] = size(ctable);
iscolmodline = cellfun(@(x) ~isempty(x), regexp(txt, 'COLOR_MODEL'));
if any(iscolmodline)
    colmodel = regexprep(txt{iscolmodline}, 'COLOR_MODEL', '');
    colmodel = strtrim(lower(regexprep(colmodel, '[#=]', '')));
else
    if nc == 8
        colmodel = 'rgb';
    elseif nc == 10
        colmodel = 'cmyk';
    else
        error('Cannot parse this format .cpt file yet');
    end
end
%     try
%         temp = str2num(strvcat(txt(~isheader & ~isfooter)));
%         if size(temp,2) == 8
%             colmodel = 'rgb';
%         elseif size(temp,2) == 10
%             colmodel = 'cmyk';
%         else % grayscale, maybe others
%             error('Cannot parse this format .cpt file yet');
%         end
%     catch % color names, mixed formats, dash placeholders
%         error('Cannot parse this format .cpt file yet');
%     end
% end
%     

% 
% iscmod = strncmp(txt, '# COLOR_MODEL', 13);
% 
% 
% if ~any(iscmod)
%     isrgb = true;
% else
%     cmodel = strtrim(regexprep(txt{iscmod}, '# COLOR_MODEL =', ''));
%     if strcmp(cmodel, 'RGB')
%         isrgb = true;
%     elseif strcmp(cmodel, 'HSV')
%         isrgb = false;
%     else
%         error('Unknown color model: %s', cmodel);
%     end
% end

% Reformat color table into one column of colors

cpt = zeros(nr*2, 4);
cpt(1:2:end,:) = ctable(:,1:4);
cpt(2:2:end,:) = ctable(:,5:8);

% Ticks

ticks = unique(cpt(:,1));

% Choose number of colors for output

if isnan(ncol)
    
    endpoints = unique(cpt(:,1));
    
    % For gradient-ed blocks, assure at least 4 steps between endpoints
    
    issolid = all(ctable(:,2:4) == ctable(:,6:8), 2);
    
    for ie = 1:length(issolid)
        if ~issolid(ie)
            temp = linspace(endpoints(ie), endpoints(ie+1), 11)';
            endpoints = [endpoints; temp(2:end-1)];
        end
    end
    
    endpoints = sort(endpoints);
    
    % Determine largest step size that resolves all endpoints
    
    space = diff(endpoints);
    space = unique(space);
    space = roundn(space, -3); % To avoid floating point issues when converting to integers
    
    nspace = length(space);
    if ~isscalar(space)
        
        fac = 1;
        tol = .001;
        while 1
            if all(space >= 1 & (space - round(space)) < tol)
                space = round(space);
                break;
            else
                space = space * 10;
                fac = fac * 10;
            end
        end
        
        pairs = nchoosek(space, 2);
        np = size(pairs,1);
        commonsp = zeros(np,1);
        for ip = 1:np
            commonsp(ip) = gcd(pairs(ip,1), pairs(ip,2));
        end
        
        space = min(commonsp);
        space = space/fac;
    end
            
    ncol = (max(endpoints) - min(endpoints))./space;
    ncol = min(ncol, 256);
    
end

% Remove replicates and mimic sharp breaks

isrep =  [false; ~any(diff(cpt),2)];
cpt = cpt(~isrep,:);

difc = diff(cpt(:,1));
minspace = min(difc(difc > 0));
isbreak = [false; difc == 0];
cpt(isbreak,1) = cpt(isbreak,1) + .01*minspace;

% Parse background, foreground, and nan colors

footer = txt(isfooter);
bfncol = nan(3,3);
for iline = 1:length(footer)
    if strcmp(footer{iline}(1), 'B')
        bfncol(1,:) = str2num(regexprep(footer{iline}, 'B', ''));
    elseif strcmp(footer{iline}(1), 'F')
        bfncol(2,:) = str2num(regexprep(footer{iline}, 'F', ''));
    elseif strcmp(footer{iline}(1), 'N')
        bfncol(3,:) = str2num(regexprep(footer{iline}, 'N', ''));
    end
end

% Convert to Matlab-format colormap and color limits

lims = [min(cpt(:,1)) max(cpt(:,1))];
endpoints = linspace(lims(1), lims(2), ncol+1);
midpoints = (endpoints(1:end-1) + endpoints(2:end))/2;

cmap = interp1(cpt(:,1), cpt(:,2:4), midpoints);
switch colmodel
    case 'rgb'
        cmap = cmap ./ 255;
        bfncol = bfncol ./ 255;
    case 'hsv'
        cmap(:,1) = cmap(:,1)./300;
        cmap = hsv2rgb(cmap);
        
        bfncol(:,1) = bfncol(:,1)./300;
        bfncol = hsv2rgb(bfncol);
        
    case 'cmyk'
        error('CMYK color conversion not yet supported');
end

%------------------------------
% Subfunction: Display all
% colormaps
%------------------------------c

function plotcmaps(folder)

Files = dir(fullfile(folder, '*.cpt'));
nfile = length(Files);
ncol = 3; %ceil(nfile/nr);
nr = ceil(nfile/ncol);%3;
width = (1 - .05*2)/ncol;
height = (1-.05*2)/nr;
left = .05 + (0:ncol-1)*width;
bot = .05 + (0:nr-1)*height;

[l, b] = meshgrid(left, bot);
w = width * .8;
h = height * .5;

figure;
for ifile = 1:nfile
    ax(ifile) = subplot('position', [l(ifile), b(ifile), width, height]);
    cptcmap(Files(ifile).name, 'mapping', 'direct');
    cb(ifile) = colorbar('horiz', 'position', [l(ifile), b(ifile), w, h]);
    set(cb(ifile), 'fontsize', 8);
    cbfreeze(cb(ifile), 'on');
    text(0, .7, Files(ifile).name, 'interpreter', 'none', 'fontsize', 10);
end

set(ax, 'ylim', [0 1], 'xlim', [0 1], 'visible', 'off');

%****** Extra for FEX

function [Param, extra] = parsepv(Param, pvpairs, varargin)
%PARSEPV Parses parameter/value pairs
%
% NewParam = parsepv(Param, pvpairs)
% [NewParam, extra] = parsepv(Param, pvpairs, 'returnextra')
%
% This function is an extension of parse_pv_pairs.  It allows the option of
% returning unrecognized parameter/value pairs, rather than erroring.
%
% Input variables:
%
%   Param:          1 x 1 structure holding default parameters (fieldnames)
%                   and values 
%
%   pvpairs:        1 x n cell array of parameter/value pairs
%
%   'returnextra':  if this string is included, the function will return a
%                   cell array holding any unrecognized parameters and the
%                   corresponding values.  Otherwise, it will error if a
%                   parameter is not recognized.  
%
% Output variables:
%
%   NewParam:       1 x 1 struct identical to Param but with defaults
%                   replaced by the values from pvpairs  
%
%   extra:          1 x m cell array of any unrecognized parameter/value
%                   pairs

% Copyright 2009 Kelly Kearney

if nargin == 3
    returnextra = strcmpi(varargin{1}, 'returnextra');
else
    returnextra = false;
end

npv = length(pvpairs);
n = npv/2;

if n~=floor(n)
  error 'Property/value pairs must come in PAIRS.'
end
if n<=0
  % just return the defaults
  if returnextra
      extra = cell(0);
  end
  return
end

if ~isstruct(Param)
    error 'No structure for defaults was supplied'
end


if returnextra
    extra = cell(0);
end

% there was at least one pv pair. process any supplied
propnames = fieldnames(Param);
lpropnames = lower(propnames);
for i=1:n
    p_i = lower(pvpairs{2*i-1});
    v_i = pvpairs{2*i};
  
    ind = strmatch(p_i,lpropnames,'exact');
    if isempty(ind)
        ind = find(strncmp(p_i,lpropnames,length(p_i)));
        if isempty(ind)
            if returnextra
                extra = [extra pvpairs(2*i-1:2*i)];
                continue;
            else
                error(['No matching property found for: ',pvpairs{2*i-1}]);
            end
        elseif length(ind)>1
            error(['Ambiguous property name: ',pvpairs{2*i-1}])
        end
    end
    p_i = propnames{ind};

    % override the corresponding default in params
    Param = setfield(Param,p_i,v_i);
  
end


